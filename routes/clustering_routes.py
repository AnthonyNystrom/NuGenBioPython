"""
Routes for clustering analysis
"""
from flask import Blueprint, request, jsonify
import base64
import math
from io import BytesIO

from dependencies import np, KMeans, DBSCAN, AgglomerativeClustering, linkage, dendrogram, plt, Cluster
from utils.plot_helpers import (
    svg_markup,
    style_axes, set_title,
    LABEL_COLOR, MUTED_COLOR, AXIS_COLOR,
    subplots as oo_subplots,
)
from utils.request_helpers import clamp_int, clamp_float, error_response

bp = Blueprint('clustering', __name__, url_prefix='/api')


@bp.route('/clustering/analyze', methods=['POST'])
def clustering_analyze():
    try:
        data = request.get_json(silent=True) or {}
        matrix = np.array(data.get('matrix', []))
        method = data.get('method', 'kmeans')
        n_clusters = int(data.get('n_clusters', 3))

        if matrix.size == 0:
            return jsonify({'success': False, 'error': 'No data matrix provided'})

        results = {}

        if method == 'kmeans':
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            labels = kmeans.fit_predict(matrix)
            results = {
                'labels': labels.tolist(),
                'centers': kmeans.cluster_centers_.tolist(),
                'inertia': float(kmeans.inertia_)
            }

        elif method == 'hierarchical':
            linkage_matrix = linkage(matrix, method='ward')

            # Pick a color threshold that carves out the requested number of
            # clusters. scipy's dendrogram uses `color_threshold` to color
            # branches above/below that merge distance.
            n = matrix.shape[0]
            if n_clusters and n_clusters > 0 and n > n_clusters:
                # Distance at which the tree would be cut to yield n_clusters
                sorted_d = sorted(linkage_matrix[:, 2])
                color_threshold = sorted_d[-(n_clusters - 1)] if n_clusters > 1 else sorted_d[-1] + 1
            else:
                color_threshold = None

            # Figure size scales with leaf count. Wider for >30 samples,
            # taller only at extreme counts (the dendrogram is naturally
            # horizontal).
            fig_w = max(9, min(22, 7 + math.sqrt(max(n, 1)) * 0.7))
            fig_h = max(5, min(10, 4 + math.log2(max(n, 2)) * 0.5))
            fig, ax = oo_subplots(figsize=(fig_w, fig_h), dpi=100)

            dendrogram(
                linkage_matrix, ax=ax,
                color_threshold=color_threshold,
                above_threshold_color=MUTED_COLOR,
                leaf_font_size=max(7, 10 - max(0, n - 20) * 0.05),
                leaf_rotation=90 if n > 15 else 0,
            )
            if color_threshold is not None:
                ax.axhline(y=color_threshold, color='#ef4444', linestyle='--',
                           linewidth=0.8, alpha=0.7, zorder=1)
                # Place label inside the axes so it isn't clipped.
                x_right = ax.get_xlim()[1]
                ax.text(x_right - (x_right - ax.get_xlim()[0]) * 0.01,
                        color_threshold,
                        f'cut → {n_clusters} clusters',
                        va='bottom', ha='right',
                        fontsize=8, color='#ef4444')

            set_title(ax, 'Hierarchical Clustering Dendrogram', pad=10)
            ax.set_xlabel('Sample', fontsize=9, color=LABEL_COLOR)
            ax.set_ylabel('Distance', fontsize=9, color=LABEL_COLOR)
            style_axes(ax)

            dendro_url = svg_markup(fig, title='Clustering dendrogram')

            clustering = AgglomerativeClustering(n_clusters=n_clusters)
            labels = clustering.fit_predict(matrix)

            results = {
                'labels': labels.tolist(),
                'linkage_matrix': linkage_matrix.tolist(),
                'dendrogram': dendro_url,
                'color_threshold': color_threshold,
            }

        elif method == 'dbscan':
            eps = clamp_float(data.get('eps'), 0.5, lo=1e-9)
            min_samples = clamp_int(data.get('min_samples'), 5, lo=1)

            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            labels = dbscan.fit_predict(matrix)

            results = {
                'labels': labels.tolist(),
                'n_clusters': len(set(labels)) - (1 if -1 in labels else 0),
                'n_noise': list(labels).count(-1)
            }

        return jsonify({'success': True, 'results': results, 'method': method})
    except Exception as e:
        return error_response(e, context='clustering_routes.clustering_analyze')


# Bio.Cluster endpoints
@bp.route('/clustering/biocluster/kcluster', methods=['POST'])
def biocluster_kmeans():
    """Bio.Cluster k-means clustering"""
    try:
        data = request.get_json(silent=True) or {}
        matrix = np.array(data.get('matrix', []), dtype=float)
        n_clusters = int(data.get('n_clusters', 3))
        n_pass = int(data.get('n_pass', 10))
        dist = data.get('distance', 'e')  # e=Euclidean, c=Pearson correlation

        if matrix.size == 0:
            return jsonify({'success': False, 'error': 'No data matrix provided'})

        # Run Bio.Cluster k-means
        clusterid, error, nfound = Cluster.kcluster(matrix, nclusters=n_clusters, npass=n_pass, dist=dist)

        # Get cluster centers
        cdata, cmask = Cluster.clustercentroids(matrix, clusterid=clusterid)

        return jsonify({
            'success': True,
            'method': 'Bio.Cluster.kcluster',
            'results': {
                'labels': clusterid.tolist(),
                'centers': cdata.tolist(),
                'error': float(error),
                'nfound': int(nfound),
                'n_clusters': n_clusters
            }
        })
    except Exception as e:
        return error_response(e, context='clustering_routes.biocluster_kmeans')


@bp.route('/clustering/biocluster/kmedoids', methods=['POST'])
def biocluster_kmedoids():
    """Bio.Cluster k-medoids clustering"""
    try:
        data = request.get_json(silent=True) or {}
        matrix = np.array(data.get('matrix', []), dtype=float)
        n_clusters = int(data.get('n_clusters', 3))
        n_pass = int(data.get('n_pass', 10))
        dist = data.get('distance', 'e')

        if matrix.size == 0:
            return jsonify({'success': False, 'error': 'No data matrix provided'})

        # Calculate distance matrix first
        distmatrix = Cluster.distancematrix(matrix, dist=dist)

        # Run k-medoids
        clusterid, error, nfound = Cluster.kmedoids(distmatrix, nclusters=n_clusters, npass=n_pass)

        return jsonify({
            'success': True,
            'method': 'Bio.Cluster.kmedoids',
            'results': {
                'labels': clusterid.tolist(),
                'error': float(error),
                'nfound': int(nfound),
                'n_clusters': n_clusters
            }
        })
    except Exception as e:
        return error_response(e, context='clustering_routes.biocluster_kmedoids')


@bp.route('/clustering/biocluster/treecluster', methods=['POST'])
def biocluster_hierarchical():
    """Bio.Cluster hierarchical clustering"""
    try:
        data = request.get_json(silent=True) or {}
        matrix = np.array(data.get('matrix', []), dtype=float)
        method = data.get('linkage', 's')  # s=single, m=complete, a=average, c=centroid
        dist = data.get('distance', 'e')  # e=Euclidean, c=correlation

        if matrix.size == 0:
            return jsonify({'success': False, 'error': 'No data matrix provided'})

        # Run hierarchical clustering
        tree = Cluster.treecluster(matrix, method=method, dist=dist)

        # Convert tree to list format
        tree_data = []
        for i in range(len(tree)):
            node = tree[i]
            tree_data.append({
                'left': int(node.left),
                'right': int(node.right),
                'distance': float(node.distance)
            })

        return jsonify({
            'success': True,
            'method': 'Bio.Cluster.treecluster',
            'results': {
                'tree': tree_data,
                'n_nodes': len(tree_data)
            }
        })
    except Exception as e:
        return error_response(e, context='clustering_routes.biocluster_hierarchical')


@bp.route('/clustering/biocluster/somcluster', methods=['POST'])
def biocluster_som():
    """Bio.Cluster Self-Organizing Map clustering"""
    try:
        data = request.get_json(silent=True) or {}
        matrix = np.array(data.get('matrix', []), dtype=float)
        xdim = int(data.get('xdim', 3))
        ydim = int(data.get('ydim', 3))
        n_iter = int(data.get('n_iter', 1000))

        if matrix.size == 0:
            return jsonify({'success': False, 'error': 'No data matrix provided'})

        # Run SOM clustering
        clusterid, celldata = Cluster.somcluster(matrix, nxgrid=xdim, nygrid=ydim, niter=n_iter)

        return jsonify({
            'success': True,
            'method': 'Bio.Cluster.somcluster',
            'results': {
                'labels': clusterid.tolist(),
                'grid_shape': [xdim, ydim],
                'celldata': celldata.tolist()
            }
        })
    except Exception as e:
        return error_response(e, context='clustering_routes.biocluster_som')


@bp.route('/clustering/biocluster/pca', methods=['POST'])
def biocluster_pca():
    """Bio.Cluster Principal Component Analysis"""
    try:
        data = request.get_json(silent=True) or {}
        matrix = np.array(data.get('matrix', []), dtype=float)

        if matrix.size == 0:
            return jsonify({'success': False, 'error': 'No data matrix provided'})

        # Run PCA
        columnmean, coordinates, components, eigenvalues = Cluster.pca(matrix)

        return jsonify({
            'success': True,
            'method': 'Bio.Cluster.pca',
            'results': {
                'coordinates': coordinates.tolist(),
                'components': components.tolist(),
                'eigenvalues': eigenvalues.tolist(),
                'columnmean': columnmean.tolist(),
                'explained_variance': [float(ev / sum(eigenvalues) * 100) for ev in eigenvalues]
            }
        })
    except Exception as e:
        return error_response(e, context='clustering_routes.biocluster_pca')


@bp.route('/clustering/biocluster/distancematrix', methods=['POST'])
def biocluster_distance():
    """Bio.Cluster distance matrix calculation"""
    try:
        data = request.get_json(silent=True) or {}
        matrix = np.array(data.get('matrix', []), dtype=float)
        dist = data.get('distance', 'e')  # e=Euclidean, c=correlation, a=absolute correlation, etc.

        if matrix.size == 0:
            return jsonify({'success': False, 'error': 'No data matrix provided'})

        # Calculate distance matrix
        distmatrix = Cluster.distancematrix(matrix, dist=dist)

        # Convert to full matrix format
        n = len(matrix)
        full_matrix = [[0.0] * n for _ in range(n)]
        for i in range(n):
            for j in range(i):
                full_matrix[i][j] = distmatrix[i][j]
                full_matrix[j][i] = distmatrix[i][j]

        return jsonify({
            'success': True,
            'method': 'Bio.Cluster.distancematrix',
            'results': {
                'distance_matrix': full_matrix,
                'distance_method': dist,
                'shape': [n, n]
            }
        })
    except Exception as e:
        return error_response(e, context='clustering_routes.biocluster_distance')
