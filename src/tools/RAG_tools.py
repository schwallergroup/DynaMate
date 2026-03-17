from paperqa import Docs
from tqdm import tqdm
import os
import os.path
import contextlib
import pickle
from dotenv import load_dotenv, set_key
from src import constants

load_dotenv(dotenv_path=constants.ENV_FILE)

documents: Docs | None = None


def _ensure_openai_key() -> None:
    key = os.environ.get("OPENAI_API_KEY")
    if not key:
        print("--- Missing API Key: OPENAI_API_KEY (required for paper search) ---")
        key = input("Please enter your OpenAI API key: ").strip()
        if key:
            os.environ["OPENAI_API_KEY"] = key
            set_key(str(constants.ENV_FILE), "OPENAI_API_KEY", key)


def _load_documents() -> Docs:
    _ensure_openai_key()
    docs = Docs()

    pdf_files = list(constants.PAPER_DIR.rglob("*.pdf")) if constants.PAPER_DIR.exists() else []
    total_files = len(pdf_files)

    if total_files == 0:
        return None

    pickled_docs = "my_docs.pkl"

    if not os.path.exists(pickled_docs):
        with tqdm(
            total=total_files,
            desc="Loading PDFs",
            bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} ({percentage:.1f}%) [ time left: {remaining}, time spent: {elapsed}]",
        ) as pbar:
            for file_path in pdf_files:
                with open(os.devnull, "w") as fnull:
                    with contextlib.redirect_stdout(fnull), contextlib.redirect_stderr(fnull):
                        docs.add(file_path)
                pbar.update(1)
        with open(pickled_docs, "wb") as f:
            pickle.dump(docs, f)
    else:
        with open(pickled_docs, "rb") as f:
            docs = pickle.load(f)

    return docs


def search_papers(query: dict):
    global documents
    _ensure_openai_key()

    if not documents:
        documents = _load_documents()

    if documents is None:
        return (f"No PDF files found in {constants.PAPER_DIR}. "
                "Paper search is unavailable. Continuing without literature context.")

    if isinstance(query, dict):
        query = query.get("query")

    if not isinstance(query, str):
        return f"search_papers expected a string query, got: {type(query)}"

    result = documents.query(query)
    answer = result.formatted_answer
    if "I cannot answer." in answer:
        answer += f" Check to ensure there's papers in {constants.PAPER_DIR}"
    return answer
