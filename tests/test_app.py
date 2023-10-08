from vcftoapi.server import create_app

def test_index(client):
    response = client.get("/")
    assert response.data == b"vcftoapi server is running"


def test_index(client):
    response = client.get("/")
    assert response.data == b"vcftoapi server is running"


def test_configuration(client):
    response = client.get("/configuration")
    assert response.data != None

def test_chromosomes(client):
    response = client.get("/chromosomes")
    assert response.data != None


def test_samples(client):
    response = client.get("/samples")
    assert response.data != None


def test_genes(client):
    response = client.get("/genes")
    assert response.data != None

def test_chromosomes(client):
    response = client.get("/chromosomes")
    assert response.data != None

def test_genomic_window_summary(client):
    response = client.get("/genomic_window_summary")
    assert response.data != None

def test_pca(client):
    response = client.get("/pca")
    assert response.data == None

def test_phylo_cluster(client):
    response = client.get("/phylo_cluster")
    assert response.data != None

def test_variant_calls(client):
    response = client.get("/variant_calls")
    assert response.data != None
