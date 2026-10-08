""" This file contains the code for calling all LLM APIs. """

import os
from functools import partial
import tiktoken
# from schema import TooLongPromptError, LLMError

# These constants were removed from the anthropic package once the legacy
# Text Completions API was retired. We keep the strings here since this repo
# only ever used them for prompt formatting / logging, not for the API call itself.
HUMAN_PROMPT = "\n\nHuman:"
AI_PROMPT = "\n\nAssistant:"
enc = tiktoken.get_encoding("cl100k_base")

try:
    from helm.common.authentication import Authentication
    from helm.common.request import Request, RequestResult
    from helm.proxy.accounts import Account
    from helm.proxy.services.remote_service import RemoteService
    # setup CRFM API
    auth = Authentication(api_key=open("crfm_api_key.txt").read().strip())
    service = RemoteService("https://crfm-models.stanford.edu")
    account: Account = service.get_account(auth)
except Exception as e:
    print(e)
    print("Could not load CRFM API key crfm_api_key.txt.")

try:   
    import anthropic
    #setup anthropic API key
    anthropic_client = anthropic.Anthropic(api_key=open("claude_api_key.txt").read().strip())
except Exception as e:
    print(e)
    print("Could not load anthropic API key claude_api_key.txt.")

try:
    import openai
    from openai import OpenAI
    organization, api_key  =  open("openai_api_key.txt").read().strip().split(":")    
    os.environ["OPENAI_API_KEY"] = api_key 
    client = OpenAI(api_key=os.environ.get("OPENAI_API_KEY"))
except Exception as e:
    print(e)
    print("Could not load OpenAI API key openai_api_key.txt.")


def log_to_file(log_file, prompt, completion, model, max_tokens_to_sample):
    """ Log the prompt and completion to a file."""
    with open(log_file, "a") as f:
        f.write("\n===================prompt=====================\n")
        f.write(f"{HUMAN_PROMPT} {prompt} {AI_PROMPT}")
        num_prompt_tokens = len(enc.encode(f"{HUMAN_PROMPT} {prompt} {AI_PROMPT}"))
        f.write(f"\n==================={model} response ({max_tokens_to_sample})=====================\n")
        f.write(completion)
        num_sample_tokens = len(enc.encode(completion))
        f.write("\n===================tokens=====================\n")
        f.write(f"Number of prompt tokens: {num_prompt_tokens}\n")
        f.write(f"Number of sampled tokens: {num_sample_tokens}\n")
        f.write("\n\n")


def complete_text_claude(prompt, stop_sequences=[HUMAN_PROMPT], model="claude-v1", max_tokens_to_sample = 2000, temperature=0.5, log_file=None, **kwargs):
    """ Call the Claude API to complete a prompt."""

    ai_prompt = AI_PROMPT
    if "ai_prompt" in kwargs is not None:
        ai_prompt = kwargs["ai_prompt"]
        del kwargs["ai_prompt"]
    # model = "claude-2"
    if model.startswith("claude-3") or model.startswith("claude-sonnet") or model.startswith("claude-opus") or model.startswith("claude-haiku"):
        messages = [
            {'role': 'user', 'content': f"{HUMAN_PROMPT} {prompt}"}
        ]
        rsp = anthropic_client.messages.create(
            model=model,
            messages=messages,
            max_tokens=max_tokens_to_sample
        )
        completion = rsp.content[0].text
        if log_file is not None:
            log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
        return completion
    try:
        rsp = anthropic_client.completions.create(
            prompt=f"{HUMAN_PROMPT} {prompt} {ai_prompt}",
            stop_sequences=stop_sequences,
            model=model,
            temperature=temperature,
            max_tokens_to_sample=max_tokens_to_sample,
            **kwargs
        )
    except anthropic.APIStatusError as e:
        print(e)
        exit()
        raise TooLongPromptError()
    except Exception as e:
        exit()
        raise LLMError(e)

    completion = rsp.completion
    if log_file is not None:
        log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
    return completion


def get_embedding_crfm(text, model="openai/gpt-4-0314"):
    request = Request(model="openai/text-similarity-ada-001", prompt=text, embedding=True)
    request_result: RequestResult = service.make_request(auth, request)
    return request_result.embedding 

def complete_text_crfm(prompt=None, stop_sequences = None, model="openai/gpt-4-0314",  max_tokens_to_sample=2000, temperature = 0.5, log_file=None, messages = None, **kwargs):

    random = log_file
    if messages:
        request = Request(
                prompt=prompt, 
                messages=messages,
                model=model, 
                stop_sequences=stop_sequences,
                temperature = temperature,
                max_tokens = max_tokens_to_sample,
                random = random
            )
    else:
        print("model", model)
        print("max_tokens", max_tokens_to_sample)
        request = Request(
                prompt=prompt, 
                model=model, 
                stop_sequences=stop_sequences,
                temperature = temperature,
                max_tokens = max_tokens_to_sample,
                random = random
        )

    try:      
        request_result: RequestResult = service.make_request(auth, request)
    except Exception as e:
        # probably too long prompt
        print(e)
        exit()
        # raise TooLongPromptError()

    if request_result.success == False:
        print(request.error)
        # raise LLMError(request.error)
    completion = request_result.completions[0].text
    if log_file is not None:
        log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
    return completion


def complete_text_openai(prompt, stop_sequences=[], model="gpt-3.5-turbo", max_tokens_to_sample=2000, temperature=0.5, log_file=None, **kwargs):

    """ Call the OpenAI API to complete a prompt."""
    raw_request = {
          "model": model,
        #   "temperature": temperature,
        #   "max_completion_tokens": max_tokens_to_sample,
        #   "stop": stop_sequences or None,  # API doesn't like empty list
          **kwargs
    }
    if model.startswith("gpt-3.5") or model.startswith("gpt-4") or model.startswith("o1"):
        # Requires openai==1.42.0
        messages = [{"role": "user", "content": prompt}]
        response = client.chat.completions.create(**{"messages": messages,**raw_request})
        completion = response.choices[0].message.content
    else:
        response = client.completions.create(**{"prompt": prompt,**raw_request})
        completion = response.choices[0].text
    if log_file is not None:
        log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
    return completion

# ---- Local model support (Ollama / any OpenAI-compatible server) ----
# ---- OpenRouter support (OpenAI-compatible endpoint, hosted, free tier) ----
try:
    from openai import OpenAI as _OpenAI_openrouter
    openrouter_client = _OpenAI_openrouter(
        base_url="https://openrouter.ai/api/v1",
        api_key=open("openrouter_api_key.txt").read().strip(),
        max_retries=0,
    )
except Exception as e:
    print(e)
    print("Could not load OpenRouter API key openrouter_api_key.txt.")
# ---- Gemini support (OpenAI-compatible endpoint, hosted, free tier) ----
try:
    from openai import OpenAI as _OpenAI_gemini
    gemini_client = _OpenAI_gemini(
    base_url="https://generativelanguage.googleapis.com/v1beta/openai/",
    api_key=open("gemini_api_key.txt").read().strip(),
    max_retries=0,
    timeout=120.0,  # fail fast rather than silently hanging for the 10-minute default
    )
except Exception as e:
    print(e)
    print("Could not load Gemini API key gemini_api_key.txt.")

def complete_text_gemini(prompt, stop_sequences=None, model="gemini-3.5-flash-lite",
                         max_tokens_to_sample=8000, temperature=0.5, log_file=None,
                         max_retries=8, reasoning_effort="low", **kwargs):
    import openai as _openai_module
    import time
    for attempt in range(max_retries):
        try:
            response = gemini_client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": prompt}],
                temperature=temperature,
                max_tokens=max_tokens_to_sample,
                stop=stop_sequences or None,
                reasoning_effort=reasoning_effort,
            )
            completion = response.choices[0].message.content
            if not completion or not completion.strip():
                print(f"Gemini returned an empty completion (attempt {attempt + 1}/{max_retries}), retrying...")
                continue
            if log_file is not None:
                log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
            return completion
        except _openai_module.RateLimitError as e:
            wait = min(60, 2 ** attempt)
            print(f"Gemini rate limit hit (attempt {attempt + 1}/{max_retries}): {e}. Waiting {wait}s...")
            time.sleep(wait)
        except Exception as e:
            wait = min(60, 2 ** attempt)
            print(f"Gemini error (attempt {attempt + 1}/{max_retries}): {type(e).__name__}: {e}. Waiting {wait}s...")
            time.sleep(wait)
    raise RuntimeError(f"Gemini call failed after {max_retries} retries.")
def complete_text_gemma(prompt, stop_sequences=None, model="gemma-4-31b-it",
                        max_tokens_to_sample=8000, temperature=0.5, log_file=None,
                        max_retries=8, **kwargs):
    import openai as _openai_module
    import time
    for attempt in range(max_retries):
        try:
            response = gemini_client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": prompt}],
                temperature=temperature,
                max_tokens=max_tokens_to_sample,
                stop=stop_sequences or None,
            )
            completion = response.choices[0].message.content
            if not completion or not completion.strip():
                print(f"Gemma returned an empty completion (attempt {attempt + 1}/{max_retries}), retrying...")
                continue
            if log_file is not None:
                log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
            return completion
        except _openai_module.RateLimitError as e:
            wait = min(60, 2 ** attempt)
            print(f"Gemma rate limit hit (attempt {attempt + 1}/{max_retries}): {e}. Waiting {wait}s...")
            time.sleep(wait)
        except Exception as e:
            # Broad catch: Gemini's OpenAI-compat layer has returned at least one
            # malformed/non-standard error (500 INTERNAL) that may not surface as
            # a normal openai.APIError. Retry regardless of the exact exception type.
            wait = min(60, 2 ** attempt)
            print(f"Gemma error (attempt {attempt + 1}/{max_retries}): {type(e).__name__}: {e}. Waiting {wait}s...")
            time.sleep(wait)
    raise RuntimeError(f"Gemma call failed after {max_retries} retries.")
def complete_text_openrouter(prompt, stop_sequences=None, model="meta-llama/llama-3.3-70b-instruct:free",
                             max_tokens_to_sample=2000, temperature=0.5, log_file=None,
                             max_retries=8, **kwargs):
    import openai as _openai_module
    import time
    for attempt in range(max_retries):
        try:
            response = openrouter_client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": prompt}],
                temperature=temperature,
                max_tokens=max_tokens_to_sample,
                stop=stop_sequences or None,
            )
            completion = response.choices[0].message.content
            if log_file is not None:
                log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
            return completion
        except _openai_module.RateLimitError as e:
            wait = min(60, 2 ** attempt)
            print(f"OpenRouter rate limit hit (attempt {attempt + 1}/{max_retries}), waiting {wait}s...")
            time.sleep(wait)
        except _openai_module.APIError as e:
            wait = min(60, 2 ** attempt)
            print(f"OpenRouter API error (attempt {attempt + 1}/{max_retries}): {e}. Waiting {wait}s...")
            time.sleep(wait)
    raise RuntimeError(f"OpenRouter call failed after {max_retries} retries.")
# ---- Groq support (OpenAI-compatible endpoint, hosted, free tier) ----
try:
    from openai import OpenAI as _OpenAI_groq
    groq_client = _OpenAI_groq(
        base_url="https://api.groq.com/openai/v1",
        api_key=open("groq_api_key.txt").read().strip(),
        max_retries=0,  # disable the SDK's own silent retries -- we handle retries explicitly below
    )
except Exception as e:
    print(e)
    print("Could not load Groq API key groq_api_key.txt.")


def complete_text_groq(prompt, stop_sequences=None, model="llama-4-scout-17b-16e-instruct",
                       max_tokens_to_sample=8000, temperature=0.5, log_file=None,
                       max_retries=8, **kwargs):
    # Call Groq's hosted API, with explicit retry-with-backoff on rate limits (429)
    # and transient server errors.
    import openai as _openai_module
    import time
    # gpt-oss models are reasoning models: internal "thinking" tokens count
    # against the same max_tokens budget as the visible answer. Without an
    # explicit low reasoning effort, the model can exhaust its whole budget
    # thinking and return an empty completion with no error.
    extra_body = {}
    if "gpt-oss" in model:
        extra_body["reasoning_effort"] = "low"
    elif "qwen3" in model.lower():
        extra_body["reasoning_effort"] = "none"
    for attempt in range(max_retries):
        try:
            response = groq_client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": prompt}],
                temperature=temperature,
                max_tokens=max_tokens_to_sample,
                stop=stop_sequences or None,
                extra_body=extra_body or None,
            )
            completion = response.choices[0].message.content
            if not completion or not completion.strip():
                print(f"Groq returned an empty completion (attempt {attempt + 1}/{max_retries}), retrying...")
                continue
            if log_file is not None:
                log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
            return completion
        except _openai_module.RateLimitError as e:
            wait = min(60, 2 ** attempt)
            print(f"Groq rate limit hit (attempt {attempt + 1}/{max_retries}), waiting {wait}s...")
            time.sleep(wait)
        except _openai_module.APIError as e:
            wait = min(60, 2 ** attempt)
            print(f"Groq API error (attempt {attempt + 1}/{max_retries}): {e}. Waiting {wait}s...")
            time.sleep(wait)
    raise RuntimeError(f"Groq call failed after {max_retries} retries.")
try:
    from openai import OpenAI as _OpenAI
    local_client = _OpenAI(
        base_url=os.environ.get("LOCAL_LLM_BASE_URL", "http://localhost:11434/v1"),
        api_key="ollama",  # required by the client but ignored by Ollama
    )
except Exception as e:
    print(e)
    print("Could not create local LLM client (pip install openai).")

def complete_text_local(prompt, stop_sequences=None, model="llama3.1:8b",
                        max_tokens_to_sample=2000, temperature=0.5, log_file=None, **kwargs):
    """ Call a locally served model through an OpenAI-compatible endpoint. """
    response = local_client.chat.completions.create(
        model=model,
        messages=[{"role": "user", "content": prompt}],
        temperature=temperature,
        max_tokens=max_tokens_to_sample,
        stop=stop_sequences or None,
    )
    completion = response.choices[0].message.content
    if log_file is not None:
        log_to_file(log_file, prompt, completion, model, max_tokens_to_sample)
    return completion
def complete_text(prompt, log_file, model, **kwargs):
    """ Complete text using the specified model with appropriate API. """
    if model.startswith("gemma:"):
        completion = complete_text_gemma(prompt, stop_sequences=["Observation:"], log_file=log_file, model=model[len("gemma:"):], **kwargs)
    elif model.startswith("gemini:"):
        completion = complete_text_gemini(prompt, stop_sequences=["Observation:"], log_file=log_file, model=model[len("gemini:"):], **kwargs)
    elif model.startswith("openrouter:"):
        completion = complete_text_openrouter(prompt, stop_sequences=["Observation:"], log_file=log_file, model=model[len("openrouter:"):], **kwargs)
    elif model.startswith("local:"):
        completion = complete_text_local(prompt, stop_sequences=["Observation:"], log_file=log_file, model=model[len("local:"):], **kwargs)
    elif model.startswith("groq:"):
        completion = complete_text_groq(prompt, stop_sequences=["Observation:"], log_file=log_file, model=model[len("groq:"):], **kwargs)
    elif model.startswith("claude"):
        # use anthropic API
        completion = complete_text_claude(prompt, stop_sequences=[HUMAN_PROMPT, "Observation:"], log_file=log_file, model=model, **kwargs)
    elif "/" in model:
        # use CRFM API since this specifies organization like "openai/..."
        completion = complete_text_crfm(prompt, stop_sequences=["Observation:"], log_file=log_file, model=model, **kwargs)
    else:
        # use OpenAI API
        completion = complete_text_openai(prompt, stop_sequences=["Observation:"], log_file=log_file, model=model, **kwargs)
    return completion

# specify fast models for summarization etc
FAST_MODEL = "claude-v1"
def complete_text_fast(prompt, **kwargs):
    return complete_text(prompt = prompt, model = FAST_MODEL, temperature =0.01, **kwargs)
