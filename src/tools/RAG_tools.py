from paperqa import Docs, Settings
from tqdm import tqdm
import os
import os.path
import contextlib
import pickle
from src import constants


documents: Docs | None = None

def _load_documents() -> Docs:
    docs = Docs()

    pdf_files = list(constants.PAPER_DIR.glob("*.pdf"))
    total_files = len(pdf_files)

    pickled_docs = "my_docs.pkl"

    if not os.path.exists(pickled_docs):
        with tqdm(
            total=total_files,
            desc="Loading PDFs",
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} ({percentage:.1f}%) [ time left: {remaining}, time spent: {elapsed}]",
        ) as pbar:
            for file_path in pdf_files:
                # suppress output from docs.add()
                with open(os.devnull, "w") as fnull:
                    with contextlib.redirect_stdout(fnull), contextlib.redirect_stderr(fnull):
                        docs.add(file_path)
                pbar.update(1)
        # save embeddings of the documents
        with open(pickled_docs, "wb") as f:
            pickle.dump(docs, f)
    else:
        # load embeddings of the documents
        with open(pickled_docs, "rb") as f:
            docs = pickle.load(f)

    return docs

def search_papers(query: dict):
    global documents 

    if not documents:
        documents = _load_documents()

    if isinstance(query, dict):
        query = query.get("query")

    if not isinstance(query, str):
        raise ValueError(f"search_papers expected a string query, got: {type(query)}")

    paper_directory = constants.PAPER_DIR
    if paper_directory is None:
        raise ValueError(
            "'paper_dir' is None. To use this tool, the user must provide a directory with PDFs at the start."
        )
    
    settings = Settings(
        # Retrieval size — more is NOT always better
        evidence_k=8,                   # retrieve top 8 chunks per query
        max_chunk_size=800,             # avoid huge chunks; MD details are often local
        rerank_k=20,                    # lightly expand initial search before reranking

        # LLM settings
        temperature=0.1,                # scientific, deterministic tone
        answer_temperature=0.0,         # final answers must be strict, non-creative

        # Trust/scientific correctness
        require_citations=True,         # every claim must have a supporting doc
        max_tokens=4096,                # MD methods can be verbose
        summary_length=5,               # keep evidence chunks tight

        # Error-handling / agent behavior
        cohere_reranker=False,          # use built-in reranker (fast, good enough)
        retries=2,                      # avoid failures during batch queries
        timeout=120,                     # MD queries can be long
    )

    result = documents.query(query, settings=settings)
    answer = result.formatted_answer
    if "I cannot answer." in answer:
        answer += f" Check to ensure there's papers in {paper_directory}"
    return answer
