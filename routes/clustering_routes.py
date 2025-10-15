"""
Routes for clustering analysis
"""
from flask import Blueprint, request, jsonify
import base64
from io import BytesIO

from dependencies import np, KMeans, DBSCAN, AgglomerativeClustering, linkage, dendrogram, plt, Cluster

bp = Blueprint('clustering', __name__, url_prefix='/api')


@bp.route('/clustering/analyze', methods=['POST'])
def clustering_analyze():
    try:
        data = request.json
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

            # Create dendrogram
            fig, ax = plt.subplots(figsize=(10, 6))
            dendrogram(linkage_matrix, ax=ax)
            ax.set_title('Hierarchical Clustering Dendrogram')

            img_buffer = BytesIO()
            plt.savefig(img_buffer, format='png', bbox_inches='tight', dpi=150)
            img_buffer.seek(0)
            dendro_base64 = base64.b64encode(img_buffer.getvalue()).decode()
            plt.close()

            # Get cluster labels
            clustering = AgglomerativeClustering(n_clusters=n_clusters)
            labels = clustering.fit_predict(matrix)

            results = {
                'labels': labels.tolist(),
                'linkage_matrix': linkage_matrix.tolist(),
                'dendrogram': f'data:image/png;base64,{dendro_base64}'
            }

        elif method == 'dbscan':
            eps = float(data.get('eps', 0.5))
            min_samples = int(data.get('min_samples', 5))

            dbscan = DBSCAN(eps=eps, min_samples=min_samples)
            labels = dbscan.fit_predict(matrix)

            results = {
                'labels': labels.tolist(),
                'n_clusters': len(set(labels)) - (1 if -1 in labels else 0),
                'n_noise': list(labels).count(-1)
            }

        return jsonify({'success': True, 'results': results, 'method': method})
    except Exception as e:
        return jsonify({'success': False, 'error': str(e)})


# Bio.Cluster endpoints
@bp.route('/clustering/biocluster/kcluster', methods=['POST'])
def biocluster_kmeans():
    """Bio.Cluster k-means clustering"""
    try:
        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/clustering/biocluster/kmedoids', methods=['POST'])
def biocluster_kmedoids():
    """Bio.Cluster k-medoids clustering"""
    try:
        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/clustering/biocluster/treecluster', methods=['POST'])
def biocluster_hierarchical():
    """Bio.Cluster hierarchical clustering"""
    try:
        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/clustering/biocluster/somcluster', methods=['POST'])
def biocluster_som():
    """Bio.Cluster Self-Organizing Map clustering"""
    try:
        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/clustering/biocluster/pca', methods=['POST'])
def biocluster_pca():
    """Bio.Cluster Principal Component Analysis"""
    try:
        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})


@bp.route('/clustering/biocluster/distancematrix', methods=['POST'])
def biocluster_distance():
    """Bio.Cluster distance matrix calculation"""
    try:
        data = request.json
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
        return jsonify({'success': False, 'error': str(e)})
