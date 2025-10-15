"""
Flask application configuration
"""
import os


def configure_app(app):
    """Configure Flask application with settings"""
    app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'nugenbiopython-secret-key-change-in-production')
    app.config['UPLOAD_FOLDER'] = 'uploads'
    app.config['MAX_CONTENT_LENGTH'] = int(os.environ.get('MAX_CONTENT_LENGTH', 16 * 1024 * 1024))  # 16MB max file size

    # Create upload folder
    os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)

    # Create example GenePop file for demonstrations
    example_file_path = os.path.join(app.config['UPLOAD_FOLDER'], 'example_popgen.gen')
    if not os.path.exists(example_file_path):
        example_content = """Example Population Genetics Study - Three Populations
Locus1
Locus2
Locus3
POP
Ind1_1, 0101 0202 0101
Ind1_2, 0102 0201 0102
Ind1_3, 0101 0202 0101
Ind1_4, 0102 0201 0102
Ind1_5, 0101 0202 0101
POP
Ind2_1, 0202 0101 0202
Ind2_2, 0201 0102 0201
Ind2_3, 0202 0101 0202
Ind2_4, 0201 0102 0201
Ind2_5, 0202 0101 0202
POP
Ind3_1, 0101 0101 0101
Ind3_2, 0102 0102 0102
Ind3_3, 0101 0101 0101
Ind3_4, 0102 0102 0102
Ind3_5, 0101 0101 0101
"""
        with open(example_file_path, 'w') as f:
            f.write(example_content)
