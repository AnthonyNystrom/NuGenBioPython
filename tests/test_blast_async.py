"""Asynchronous BLAST submission.

A BLAST against NCBI routinely runs 30-180 seconds. Running it inline held a
worker open for the whole search and put the result at the mercy of every
proxy and browser timeout in between — if any of them gave up first the user
got nothing, even though the search had completed.

/blast/submit returns a job id immediately; /blast/job/<id> is polled.
"""
import time

import pytest

import routes.blast_routes as blast_routes


VALID_XML = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput><BlastOutput_program>blastn</BlastOutput_program>
<BlastOutput_version>BLASTN 2.0</BlastOutput_version><BlastOutput_reference>x</BlastOutput_reference>
<BlastOutput_db>nt</BlastOutput_db><BlastOutput_query-ID>Q1</BlastOutput_query-ID>
<BlastOutput_query-def>query</BlastOutput_query-def><BlastOutput_query-len>40</BlastOutput_query-len>
<BlastOutput_param><Parameters><Parameters_expect>10</Parameters_expect>
<Parameters_sc-match>1</Parameters_sc-match><Parameters_sc-mismatch>-2</Parameters_sc-mismatch>
<Parameters_gap-open>5</Parameters_gap-open><Parameters_gap-extend>2</Parameters_gap-extend>
<Parameters_filter>L</Parameters_filter></Parameters></BlastOutput_param>
<BlastOutput_iterations><Iteration><Iteration_iter-num>1</Iteration_iter-num>
<Iteration_query-ID>Q1</Iteration_query-ID><Iteration_query-def>query</Iteration_query-def>
<Iteration_query-len>40</Iteration_query-len><Iteration_hits><Hit><Hit_num>1</Hit_num>
<Hit_id>gi|1</Hit_id><Hit_def>Synthetic construct</Hit_def><Hit_accession>XM_1</Hit_accession>
<Hit_len>100</Hit_len><Hit_hsps><Hsp><Hsp_num>1</Hsp_num><Hsp_bit-score>80.0</Hsp_bit-score>
<Hsp_score>40</Hsp_score><Hsp_evalue>1e-15</Hsp_evalue><Hsp_query-from>1</Hsp_query-from>
<Hsp_query-to>40</Hsp_query-to><Hsp_hit-from>1</Hsp_hit-from><Hsp_hit-to>40</Hsp_hit-to>
<Hsp_identity>40</Hsp_identity><Hsp_positive>40</Hsp_positive><Hsp_gaps>0</Hsp_gaps>
<Hsp_align-len>40</Hsp_align-len><Hsp_qseq>ATGCGTACGTATGCGTACGTATGCGTACGTATGCGTACGT</Hsp_qseq>
<Hsp_hseq>ATGCGTACGTATGCGTACGTATGCGTACGTATGCGTACGT</Hsp_hseq>
<Hsp_midline>||||||||||||||||||||||||||||||||||||||||</Hsp_midline></Hsp></Hit_hsps></Hit>
</Iteration_hits><Iteration_stat><Statistics><Statistics_db-num>1</Statistics_db-num>
<Statistics_db-len>100</Statistics_db-len><Statistics_hsp-len>0</Statistics_hsp-len>
<Statistics_eff-space>1000</Statistics_eff-space><Statistics_kappa>0.7</Statistics_kappa>
<Statistics_lambda>1.3</Statistics_lambda><Statistics_entropy>1.1</Statistics_entropy>
</Statistics></Iteration_stat></Iteration></BlastOutput_iterations></BlastOutput>"""

SEQ = "ATGCGTACGT" * 4
FORM = {"sequence": SEQ, "program": "blastn", "database": "nt"}


@pytest.fixture
def slow_blast(monkeypatch):
    """Stand in for a real NCBI round trip, without touching the network."""
    def fake(sequence, program, database, **kwargs):
        time.sleep(0.4)
        return VALID_XML
    monkeypatch.setattr(blast_routes, "run_ncbi_blast", fake)
    return fake


def _poll(client, job_id, timeout=15):
    deadline = time.time() + timeout
    seen = []
    while time.time() < deadline:
        job = client.get(f"/api/blast/job/{job_id}").get_json()
        if job.get("status") not in seen:
            seen.append(job.get("status"))
        if job.get("status") in ("done", "error"):
            return job, seen
        time.sleep(0.1)
    pytest.fail(f"job did not finish; states seen: {seen}")


def test_submit_returns_immediately(client, slow_blast):
    """The point of the change: the search must not block the response."""
    started = time.time()
    data = client.post("/api/blast/submit", data=FORM).get_json()
    elapsed = time.time() - started
    assert data["success"] is True
    assert data["job_id"]
    assert data["status"] == "queued"
    assert elapsed < 0.3, (
        f"submit took {elapsed:.2f}s; the search itself takes 0.4s, so this "
        "is still blocking")


def test_job_runs_to_completion(client, slow_blast):
    job_id = client.post("/api/blast/submit", data=FORM).get_json()["job_id"]
    job, seen = _poll(client, job_id)
    assert job["status"] == "done", job.get("error")
    assert job["num_hits"] == 1
    assert job["results"][0]["title"].startswith("gi|1")
    assert "statistics" in job


def test_job_reports_progress_states(client, slow_blast):
    job_id = client.post("/api/blast/submit", data=FORM).get_json()["job_id"]
    _job, seen = _poll(client, job_id)
    assert "done" in seen
    assert any(s in seen for s in ("queued", "running"))


def test_a_failing_search_is_reported_not_swallowed(client, monkeypatch):
    def boom(sequence, program, database, **kwargs):
        raise RuntimeError("NCBI refused the connection")
    monkeypatch.setattr(blast_routes, "run_ncbi_blast", boom)

    job_id = client.post("/api/blast/submit", data=FORM).get_json()["job_id"]
    job, _seen = _poll(client, job_id)
    assert job["status"] == "error"
    assert job["error"]


def test_worker_errors_do_not_leak_internals(client, monkeypatch):
    """The job payload goes to the browser, so it gets the same treatment."""
    def boom(sequence, program, database, **kwargs):
        raise OSError("[Errno 2] No such file: '/srv/app/uploads/x.tmp'")
    monkeypatch.setattr(blast_routes, "run_ncbi_blast", boom)

    job_id = client.post("/api/blast/submit", data=FORM).get_json()["job_id"]
    job, _seen = _poll(client, job_id)
    assert job["status"] == "error"
    assert "uploads" not in job["error"]


def test_unknown_job_is_a_404(client):
    resp = client.get("/api/blast/job/not-a-real-job")
    assert resp.status_code == 404
    assert resp.get_json()["success"] is False


def test_submit_validates_before_queueing(client):
    """A bad request must fail fast, not create a doomed job."""
    data = client.post("/api/blast/submit",
                       data={"sequence": "", "program": "blastn"}).get_json()
    assert data["success"] is False
    assert "job_id" not in data


def test_submit_rejects_an_invalid_program(client):
    data = client.post("/api/blast/submit",
                       data={"sequence": SEQ, "program": "nonsense"}).get_json()
    assert data["success"] is False


def test_result_id_is_allocated_in_the_request_context(client, slow_blast):
    """/blast/alignment and /blast/export read session['blast_result_id'];
    the worker has no request context, so it must be set at submit time."""
    with client.session_transaction() as sess:
        sess.clear()
    client.post("/api/blast/submit", data=FORM)
    with client.session_transaction() as sess:
        assert sess.get("blast_result_id")


def test_the_blocking_endpoint_still_works(client, slow_blast):
    """Kept for callers that want a single round trip."""
    data = client.post("/api/blast/search", data=FORM).get_json()
    assert data["success"] is True
    assert data["num_hits"] == 1


def test_both_paths_accept_the_same_parameters(client, slow_blast):
    """Shared parsing means sync and async cannot drift apart."""
    sync = client.post("/api/blast/search", data=FORM).get_json()
    job_id = client.post("/api/blast/submit", data=FORM).get_json()["job_id"]
    async_job, _ = _poll(client, job_id)
    assert sync["num_hits"] == async_job["num_hits"]
    assert sync["results"][0]["title"] == async_job["results"][0]["title"]
