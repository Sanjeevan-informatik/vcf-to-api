import pytest
from vcftoapi.server import create_app

@pytest.fixture
def client():
    app = create_app()

    yield app.test_client()