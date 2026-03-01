# AI-API-skill.md
> **用途**：給 AI coding / 轉寫（prompt-to-code / code-to-code / agentic coding）時，快速查表與落地整合  
> **涵蓋**：OpenAI、Anthropic (Claude)、DeepSeek（OpenAI 相容介面）  
> **更新日期**：2026-02-27  
> **來源整併**：你上傳的《AI 雲平台 API 參考手冊》（節錄與整編）+ 前述 OpenAI/Anthropic 討論 + DeepSeek 官方文件

---

## 0. 你在寫「AI API Wrapper」時一定要先釐清的 8 件事

1. **Base URL**：不同供應商不同；有些提供 OpenAI 相容層（DeepSeek）。
2. **Auth Header**：Bearer / x-api-key / version header（Anthropic 必填）。
3. **主端點**：OpenAI 推 `Responses API`；Anthropic 用 `Messages API`；DeepSeek 用 OpenAI 相容的 `chat/completions`。
4. **Model ID vs Model 版本**：  
   - OpenAI/Anthropic/DeepSeek 都以 **model 名稱**作為最重要的選擇鍵。  
   - DeepSeek 特別提醒：`/v1` 只是相容路徑，**與模型版本無關**。
5. **上下文（context window）與最大輸出**：直接決定你能放多少程式碼/文件與單次回覆可多長。
6. **Streaming**：SSE / chunk delta 格式不完全相同，務必做抽象層。
7. **工具/函數呼叫（tool calling / function calling）**：不同模型/模式支援度不同（DeepSeek 推理模式有禁用項）。
8. **錯誤處理 + 退避重試**：401/429/5xx 需要不同策略；要記錄 request-id（若有）與 token 用量。

---

## 1) OpenAI

### 1.1 Base URL & 認證
- **Base URL**：`https://api.openai.com/v1`
- **Auth**：`Authorization: Bearer $OPENAI_API_KEY`
-（可選）Project/Org header：視你的帳號/專案設定而定

### 1.2 核心端點（建議新專案優先）
- **Responses API**：`POST /v1/responses`（統一生成介面：文字/圖像輸入、文字輸出、工具、結構化輸出等）
-（相容）Chat Completions：`POST /v1/chat/completions`（較舊介面，但大量 SDK/範例仍在使用）

> **實務建議**：如果你要做「AI coding agent（會跑工具/多步推理/可追溯輸出）」→ 優先用 Responses API。

### 1.3 常用模型族群（以「coding / agent」角度）
> 這裡列的是「你最常會用到」的族群；實際可用清單請以 `GET /v1/models` 或官方 models 頁為準。

| 任務 | 建議模型（例） | 能力重點 |
|---|---|---|
| 旗艦寫碼/代理 | `gpt-4.1` / `gpt-4.1-mini` | 最新一代；coding + agentic tasks |
| 多模態/視覺 | `gpt-4o` / `gpt-4o-mini` | 文字+圖像輸入；速度快 |
| 推理（Reasoning） | `o3` / `o4-mini` / `o1-pro` | 深度推理；可調 reasoning effort |
| 預覽版 | `gpt-4.5-preview` | 實驗性前沿模型 |
| 影像生成 | `gpt-image-1`（或後繼版本） | 文生圖/圖編輯（依模型能力） |
| 語音 STT | `gpt-4o-mini-transcribe`（或同系） | 高精度轉錄 |
| 語音 TTS | `gpt-4o-mini-tts`（或同系） | 文字轉語音、多風格 |

> 你上傳的手冊也包含 OpenAI 模型與音訊/影像範例（可當落地參考）。

### 1.4 最小可用範例：Responses API
```bash
curl https://api.openai.com/v1/responses \
  -H "Authorization: Bearer $OPENAI_API_KEY" \
  -H "Content-Type: application/json" \
  -d '{
    "model": "gpt-4.1",
    "input": "Write a clean TypeScript function to parse a .env file into an object."
  }'
```

### 1.5 「寫 Wrapper」的介面抽象建議（OpenAI）
- `model`: string
- `input`: string | list[content parts]
- `stream`: bool
- `tools`: list[tool schema]（若你要做 agent）
- `response_format`: JSON schema / strict JSON（若你要強約束輸出）
- `metadata`: trace id / user id（方便你做 observability）

---

## 2) Anthropic (Claude)

### 2.1 Base URL & 認證（必填 header）
- **Base URL**：`https://api.anthropic.com`
- **端點**：`POST /v1/messages`
- **Headers**：
  - `x-api-key: $ANTHROPIC_API_KEY`
  - `anthropic-version: 2023-06-01`
  - `content-type: application/json`

### 2.2 模型（API IDs）
> 以官方模型概覽頁為準；常見（2026 初）主力如下：

| 層級 | Claude API ID（例） | 建議用途 |
|---|---|---|
| Opus | `claude-opus-4-6` | 最複雜的 agent/coding/推理 |
| Sonnet | `claude-sonnet-4-6` / `claude-sonnet-4-5`（alias） | 速度/能力平衡，日常寫碼很常用 |
| Haiku | `claude-haiku-4-5`（alias） | 低成本/高吞吐 |

> Claude 目前模型支援文字與圖片輸入、文字輸出（Vision）；是否提供 1M context（beta）取決於你使用的 header/方案。

### 2.3 最小可用範例：Messages API
```bash
curl https://api.anthropic.com/v1/messages \
  -H "x-api-key: $ANTHROPIC_API_KEY" \
  -H "anthropic-version: 2023-06-01" \
  -H "content-type: application/json" \
  -d '{
    "model": "claude-sonnet-4-6",
    "max_tokens": 1024,
    "messages": [
      {"role":"user","content":"Refactor this Python function for readability and add type hints."}
    ]
  }'
```

### 2.4 Claude 的「思考/推理」與你寫碼的關聯
- Claude 支援延伸思考（extended thinking）/ 思考預算（budget tokens）等能力，適合：
  - 大型重構（跨檔案）
  - 安全性檢查（依賴掃描、注入風險）
  - 需要多步推理的設計題

> 實作上：把「思考設定」視為 provider-specific config（不要硬塞進通用 schema）。

---

## 3) DeepSeek（OpenAI 相容 API）

### 3.1 Base URL & 認證
DeepSeek 官方明確指出：DeepSeek API **使用與 OpenAI 兼容的 API 格式**。  
- **Base URL**：`https://api.deepseek.com`  
-（相容寫法）也可設 `base_url = https://api.deepseek.com/v1`，**但這個 `/v1` 與模型版本無關**。  
- **Auth**：`Authorization: Bearer $DEEPSEEK_API_KEY`

### 3.2 核心端點（OpenAI 風格）
- Chat Completions：`POST https://api.deepseek.com/chat/completions`

### 3.3 模型與能力（DeepSeek-V3.2）
DeepSeek 文件指出：  
- `deepseek-chat` 與 `deepseek-reasoner` 對應同一代 **DeepSeek-V3.2**，上下文 **128K**。  
- `deepseek-chat`：非思考模式  
- `deepseek-reasoner`：思考模式（會輸出 reasoning chain）

| 模型 | 模式 | 上下文 | 輸出長度 | 特色 |
|---|---|---:|---:|---|
| `deepseek-chat` | 非思考 | 128K | 預設 4K、最大 8K | 通用對話/寫碼 |
| `deepseek-reasoner` | 思考 | 128K | 預設 32K、最大 64K | 會輸出 `reasoning_content` |

### 3.4 推理模型（deepseek-reasoner）注意事項（非常重要）
官方說明（摘要）：
- 會在回覆中提供 **`reasoning_content`**（思考鏈）與 **`content`**（最終回答）。
- **下一輪 messages 不可把 `reasoning_content` 拼回去**；若你把它放進 messages，會得到 **400** 錯誤。
- 不支援：Function Calling、FIM（且多個 sampling/penalty 參數不生效或會報錯）

> **落地建議**：你的 wrapper 應把 `reasoning_content` 當成「可選的 side-channel」，不要混進一般對話上下文。

### 3.5 最小可用範例（curl）
```bash
curl https://api.deepseek.com/chat/completions \
  -H "Content-Type: application/json" \
  -H "Authorization: Bearer $DEEPSEEK_API_KEY" \
  -d '{
    "model": "deepseek-chat",
    "messages": [
      {"role":"system","content":"You are a helpful assistant."},
      {"role":"user","content":"Write a Go function to validate IPv4 addresses."}
    ],
    "stream": false
  }'
```

### 3.6 服務品質/限速/錯誤（摘要）
- 官方表示：**不限制並發量**，但高流量時會維持連線並回傳 keep-alive（非流式：空行；流式：SSE 註解 `: keep-alive`）。  
- 若 **10 分鐘**後仍未開始推理，伺服器會關閉連線。  
- 常見錯誤碼：400/401/402/422/429/500/503（含原因與建議處理）

---

## 4) 跨平台能力對照（AI-QMS 支援的主要 Provider）

> 以「AI coding / agent」常用功能來比較；實際是否可用仍需以你帳號/模型支援為準。

### 4.1 Direct API Provider

| 能力 | OpenAI | Anthropic (Claude) | Google Gemini | DeepSeek | xAI (Grok) |
|---|---|---|---|---|---|
| 主要生成端點 | `/v1/responses`（推薦） | `/v1/messages` | Gemini API | `/chat/completions`（OpenAI 相容） | `/v1/chat/completions` |
| Vision（看圖） | ✅ | ✅ | ✅ | ⚠️（chat 模式限定） | ✅ |
| Tool / Function calling | ✅ | ✅ | ✅ | ✅（chat）；❌（reasoner） | ✅ |
| 長上下文 | ✅（看模型） | ✅（200K，部分 1M beta） | ✅（1M+） | ✅（128K） | ✅（128K+） |
| 推理鏈輸出 | o-series 支援 | extended thinking | Gemini 2.5 支援 | ✅（`reasoning_content`） | 視模型 |
| LiteLLM 前綴 | `openai/` | `anthropic/` | `gemini/` | `deepseek/` | `xai/` |

### 4.2 Gateway / Router Platform

| 平台 | 特色 | 模型數 | LiteLLM 前綴 |
|---|---|---|---|
| OpenRouter | 200+ 模型統一 API | 200+ | `openrouter/` |
| Together AI | 開源模型最佳化推理 | 200+ | `together_ai/` |
| Groq | LPU 超高速推理 | 20+ | `groq/` |
| Fireworks AI | 最佳化推理引擎 | 50+ | `fireworks_ai/` |
| Deep Infra | 高性價比推理 | 50+ | `deepinfra/` |

### 4.3 Local Provider

| 平台 | 特色 | LiteLLM 前綴 |
|---|---|---|
| Ollama | 100+ 本地模型，免費 | `ollama/` |
| LM Studio | OpenAI 相容本地 API | `openai/`（改 base_url） |

---

## 5) 統一 Wrapper 的「建議資料結構」

### 5.0 推薦方式：使用 LiteLLM 統一抽象層

> **AI-QMS 專案已採用此方式**（見 `src/llm_providers.py`）。LiteLLM 提供 100+ 模型的統一 API，免自建 Adapter。

```python
# pip install litellm
from litellm import completion

# 所有 provider 用同一個 completion() 函式，只差 model 前綴
response = completion(
    model="openai/gpt-4.1",          # 或 "anthropic/claude-sonnet-4-6"
    messages=[{"role": "user", "content": "Hello"}],
    num_retries=2,                    # 內建指數退避（429/5xx）
)
```

**LiteLLM model 前綴對照**：
| Provider | 前綴格式 | 範例 |
|---|---|---|
| OpenAI | `openai/{model}` | `openai/gpt-4.1` |
| Anthropic | `anthropic/{model}` | `anthropic/claude-sonnet-4-6` |
| Google | `gemini/{model}` | `gemini/gemini-2.5-pro` |
| DeepSeek | `deepseek/{model}` | `deepseek/deepseek-chat` |
| xAI | `xai/{model}` | `xai/grok-4-0709` |
| OpenRouter | `openrouter/{org}/{model}` | `openrouter/google/gemini-3-pro-preview` |
| Together AI | `together_ai/{org}/{model}` | `together_ai/meta-llama/Llama-3.3-70B-Instruct-Turbo` |
| Groq | `groq/{model}` | `groq/llama-3.3-70b-versatile` |
| Ollama (local) | `ollama/{model}` | `ollama/qwen2.5:7b` |

### 5.1 內部通用請求（你自己的 schema）
```json
{
  "provider": "openai | anthropic | google | deepseek | xai | mistral | cohere | perplexity | openrouter | together | groq | fireworks | deepinfra | ollama | lmstudio",
  "model": "string",
  "messages": [
    {"role":"system|user|assistant","content":"..."}
  ],
  "input": "optional (OpenAI responses style)",
  "stream": false,
  "max_tokens": 2048,
  "temperature": 0.2,
  "tools": [],
  "metadata": {"trace_id":"...", "user_id":"..."}
}
```

### 5.2 Provider Adapter 規則（建議）
- **若使用 LiteLLM**（推薦）：LiteLLM 內部已處理各家 API 差異，你只需管理 model 前綴與 API key
- **若自建 Adapter**：
  - **OpenAIAdapter**
    - 優先走 `/v1/responses`；若你要相容舊 SDK，再提供 `/v1/chat/completions` fallback
  - **AnthropicAdapter**
    - 固定 `/v1/messages`；強制插入 `anthropic-version` header
  - **DeepSeekAdapter**
    - 直接用 OpenAI SDK 的 chat.completions（改 base_url）
    - 若模型是 `deepseek-reasoner`：把 `reasoning_content` 從回應拆出，且禁止回填進 messages

### 5.3 Fallback Chain 設計（建議）
當主要 provider 失敗時，依優先級嘗試備援：
1. **Direct API**（最高可靠性）：OpenAI → Anthropic → Google → DeepSeek → xAI → Mistral → Cohere → Perplexity
2. **Gateway**（多模型存取）：OpenRouter → Together → Groq → Fireworks → DeepInfra
3. **Local**（零成本）：Ollama → LM Studio

---

## 6) 可靠性工程（建議你直接照抄到 production checklist）

1. **Key 管理**：只放 server-side；用 KMS/Secret Manager；避免前端曝光。
2. **Timeout**：
   - DeepSeek 高流量會 keep-alive，別用過短 read timeout；同時要處理空行/註解。
3. **429 退避**：指數退避 + jitter；保留重試上限與可觀測紀錄。
4. **Idempotency（若有）**：對「會觸發工具/外部副作用」的請求設計去重 key。
5. **觀測性**：log：provider、model、latency、input/output tokens、stop reason、request id（若回傳）。
6. **回覆格式保證**：要求 JSON 時用「嚴格 JSON」或 schema（若供應商支援），並做 parser guard。

---

## 7) 參考（官方文件）
- OpenAI API Reference（overview / responses）
- Anthropic Claude API（overview / models）
- DeepSeek API Docs（快速開始 / 模型與價格 / 限速 / 錯誤碼 / 推理模型）

---

## 附錄 A：AI-QMS 專案的實際整合狀態

AI-QMS 已透過 LiteLLM 整合 16 個 provider（見 `src/llm_providers.py`）：
- **Direct API**（8 個）：OpenAI、Anthropic、Google Gemini、DeepSeek、xAI Grok、Mistral、Cohere、Perplexity
- **Gateway**（6 個）：OpenRouter、Requesty、Together AI、Groq、Fireworks AI、Deep Infra
- **Local**（2 個）：Ollama、LM Studio

### 已實作的關鍵特性
- ✅ LiteLLM 統一抽象層（免自建 Adapter）
- ✅ 自動 Fallback Chain（16 個 provider 依序嘗試）
- ✅ Streaming 支援（Chainlit + Gradio token-by-token）
- ✅ Vision / PDF OCR（MarkItDown → Vision LLM 三層 fallback）
- ✅ `reasoning_content` side-channel 處理（DeepSeek reasoner）
- ✅ `num_retries` 指數退避重試（429/5xx）
- ✅ Arize Phoenix LLM 觀測性追蹤
- ✅ 動態模型清單更新與快取（`model_cache.json`）

### 待擴展的供應商
你上傳的《AI 雲平台 API 參考手冊》還包含 ElevenLabs、WaveSpeed、fal.ai 等平台。
若需要加入，建議在 `DEFAULT_PROVIDERS` 中新增 provider config，LiteLLM 會自動處理 API 差異。
