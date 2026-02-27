# AI 雲平台 API 參考手冊

> **用途**：供 Claude Code 外接 AI API 服務時查閱  
> **更新日期**：2026-02-20  

---

## 目錄

1. [OpenAI](#1-openai)
2. [Anthropic (Claude)](#2-anthropic-claude)
3. [Google (Gemini)](#3-google-gemini)
4. [xAI (Grok)](#4-xai-grok)
5. [ElevenLabs (語音專家)](#5-elevenlabs)
6. [WaveSpeed AI (多模態加速聚合)](#6-wavespeed-ai)
7. [fal.ai (生成式媒體 API)](#7-falai)
8. [Together AI (開源模型推論)](#8-together-ai)
9. [跨平台功能對照表](#9-跨平台功能對照表)

---

## 1. OpenAI

| 項目 | 說明 |
|------|------|
| **API Base URL** | `https://api.openai.com/v1` |
| **認證方式** | `Authorization: Bearer sk-...` |
| **SDK** | `pip install openai` / `npm install openai` |
| **API 風格** | REST + Responses API (新) / Chat Completions (傳統) |

### 1.1 LLM 模型

| 模型 ID | 類型 | 上下文窗口 | 備註 |
|---------|------|-----------|------|
| `gpt-5.2` | 旗艦 | 256K | 最新最強，推薦預設使用 |
| `gpt-5.1` | 上一代旗艦 | 256K | 穩定版 |
| `gpt-5` | 推理模型 | 256K | 可調 reasoning effort |
| `gpt-4.1` | 高 CP 值 | 1M | 長上下文、高性價比 |
| `gpt-4.1-mini` | 輕量 | 1M | 速度快、成本低 |
| `gpt-4.1-nano` | 超輕量 | 1M | 最便宜 |
| `gpt-5-codex` | 程式碼 | — | 專為 Codex CLI 優化 |

**基本呼叫方式 (Chat Completions)**：
```python
from openai import OpenAI
client = OpenAI(api_key="sk-...")

response = client.chat.completions.create(
    model="gpt-5.2",
    messages=[
        {"role": "system", "content": "You are a helpful assistant."},
        {"role": "user", "content": "Hello!"}
    ]
)
print(response.choices[0].message.content)
```

**Responses API (新版，支援 tools)**：
```python
response = client.responses.create(
    model="gpt-5.2",
    input="Explain quantum computing",
    tools=[{"type": "web_search_preview"}]
)
```

### 1.2 文生圖

| 模型 ID | 備註 |
|---------|------|
| `gpt-image-1.5` | 最新旗艦，速度快 4x，文字渲染佳 |
| `gpt-image-1` | 穩定版 |
| `gpt-image-1-mini` | 低成本 (便宜 80%) |

```python
response = client.images.generate(
    model="gpt-image-1.5",
    prompt="A futuristic city skyline at sunset",
    size="1024x1024",
    quality="high"       # low / medium / high
)
image_url = response.data[0].url
```

### 1.3 TTS (文字轉語音)

| 模型 ID | 備註 |
|---------|------|
| `gpt-4o-mini-tts` | 可控表現力、多風格 |

```python
response = client.audio.speech.create(
    model="gpt-4o-mini-tts",
    voice="alloy",          # alloy / echo / fable / onyx / nova / shimmer
    input="Hello, welcome to the future."
)
response.stream_to_file("output.mp3")
```

### 1.4 語音轉文字 (STT)

| 模型 ID | 備註 |
|---------|------|
| `gpt-4o-mini-transcribe` | 推薦，精度最高 |
| `gpt-4o-transcribe` | 替代選項 |
| `whisper-1` | 經典款 |

```python
audio_file = open("recording.mp3", "rb")
transcript = client.audio.transcriptions.create(
    model="gpt-4o-mini-transcribe",
    file=audio_file
)
print(transcript.text)
```

---

## 2. Anthropic (Claude)

| 項目 | 說明 |
|------|------|
| **API Base URL** | `https://api.anthropic.com/v1` |
| **認證方式** | `x-api-key: sk-ant-...` + `anthropic-version: 2023-06-01` |
| **SDK** | `pip install anthropic` / `npm install @anthropic-ai/sdk` |
| **API 風格** | Messages API |
| **也可透過** | AWS Bedrock / Google Vertex AI |

### 2.1 LLM 模型

| 模型 ID | 類型 | 上下文 | 定價 (Input/Output per MTok) |
|---------|------|-------|----------------------------|
| `claude-opus-4-6` | 最強旗艦 | 200K (1M beta) | $5 / $25 |
| `claude-sonnet-4-6` | 新一代平衡 | 200K (1M beta) | $3 / $15 |
| `claude-sonnet-4-5-20250929` | 前代旗艦 | 200K (1M beta) | $3 / $15 |
| `claude-haiku-4-5-20251001` | 快速低成本 | 200K | $1 / $5 |

> **備註**：Claude 目前不提供原生文生圖、TTS 或 STT 功能。如需這些功能，請搭配其他平台 API。

**基本呼叫方式**：
```python
import anthropic
client = anthropic.Anthropic(api_key="sk-ant-...")

message = client.messages.create(
    model="claude-sonnet-4-5-20250929",
    max_tokens=1024,
    messages=[
        {"role": "user", "content": "Explain quantum computing in simple terms."}
    ]
)
print(message.content[0].text)
```

**附帶圖片 (Vision)**：
```python
import base64

with open("image.png", "rb") as f:
    image_data = base64.standard_b64encode(f.read()).decode("utf-8")

message = client.messages.create(
    model="claude-sonnet-4-5-20250929",
    max_tokens=1024,
    messages=[{
        "role": "user",
        "content": [
            {"type": "image", "source": {"type": "base64", "media_type": "image/png", "data": image_data}},
            {"type": "text", "text": "What's in this image?"}
        ]
    }]
)
```

**Extended Thinking (深度推理)**：
```python
message = client.messages.create(
    model="claude-opus-4-6",
    max_tokens=16000,
    thinking={
        "type": "enabled",
        "budget_tokens": 10000
    },
    messages=[{"role": "user", "content": "Solve this complex math problem..."}]
)
```

---

## 3. Google (Gemini)

| 項目 | 說明 |
|------|------|
| **API Base URL (AI Studio)** | `https://generativelanguage.googleapis.com/v1beta` |
| **認證方式** | `?key=API_KEY` 或 OAuth |
| **SDK** | `pip install google-genai` / `npm install @google/genai` |
| **企業版** | Google Cloud Vertex AI |

### 3.1 LLM 模型

| 模型 ID | 類型 | 上下文 | 備註 |
|---------|------|-------|------|
| `gemini-3.1-pro-preview` | 最新旗艦 (Preview) | 1M | 2026/02/19 發布 |
| `gemini-3-pro-preview` | Gemini 3 系列 | 1M | 強推理與程式碼 |
| `gemini-3-flash-preview` | 快速版 | 1M | Pro 級智能、Flash 速度 |
| `gemini-2.5-pro` | 穩定旗艦 | 1M | 正式版 |
| `gemini-2.5-flash` | 高性價比 | 1M | 快速穩定 |
| `gemini-2.5-flash-lite` | 超輕量 | 1M | 最低成本 |

**基本呼叫方式**：
```python
from google import genai

client = genai.Client(api_key="YOUR_API_KEY")

response = client.models.generate_content(
    model="gemini-2.5-flash",
    contents="Explain quantum computing"
)
print(response.text)
```

**Gemini 3 系列新功能 — thinking_level**：
```python
from google.genai import types

response = client.models.generate_content(
    model="gemini-3-pro-preview",
    contents="Complex reasoning task...",
    config=types.GenerateContentConfig(
        thinking_config=types.ThinkingConfig(thinking_level="high")  # minimal / low / medium / high
    )
)
```

### 3.2 文生圖

| 模型 ID | 備註 |
|---------|------|
| `imagen-4-ultra` | 最高品質 |
| `imagen-4` | 標準版 (GA) |
| `imagen-4-fast` | 快速版 |
| `gemini-3-pro-image-preview` | Nano Banana Pro，原生整合 |
| `gemini-2.5-flash-image` | Gemini 原生圖像生成 |

```python
from google import genai
from google.genai import types

client = genai.Client(api_key="YOUR_API_KEY")

# Imagen 4
response = client.models.generate_images(
    model="imagen-4",
    prompt="A serene Japanese garden in autumn",
    config=types.GenerateImagesConfig(
        number_of_images=1
    )
)

# Gemini 原生圖片生成
response = client.models.generate_content(
    model="gemini-2.5-flash-image",
    contents="Generate an image of a cute robot",
    config=types.GenerateContentConfig(
        response_modalities=["TEXT", "IMAGE"]
    )
)
```

### 3.3 TTS (文字轉語音)

| 模型 ID | 備註 |
|---------|------|
| `gemini-2.5-flash-preview-tts` | 低延遲，適合即時應用 |
| `gemini-2.5-pro-preview-tts` | 高品質，適合有聲書/Podcast |

```python
response = client.models.generate_content(
    model="gemini-2.5-flash-preview-tts",
    contents="Hello, this is a text to speech demo.",
    config=types.GenerateContentConfig(
        response_modalities=["AUDIO"],
        speech_config=types.SpeechConfig(
            voice_config=types.VoiceConfig(
                prebuilt_voice_config=types.PrebuiltVoiceConfig(
                    voice_name="Kore"   # 可用: Puck, Charon, Kore, Fenrir, Aoede 等
                )
            )
        )
    )
)
# 回應中 response.candidates[0].content.parts[0].inline_data 包含音頻資料
```

**多語者 TTS**：
```python
config=types.GenerateContentConfig(
    response_modalities=["AUDIO"],
    speech_config=types.SpeechConfig(
        multi_speaker_voice_config=types.MultiSpeakerVoiceConfig(
            speaker_voice_configs=[
                types.SpeakerVoiceConfig(speaker="Host", voice_config=types.VoiceConfig(
                    prebuilt_voice_config=types.PrebuiltVoiceConfig(voice_name="Kore"))),
                types.SpeakerVoiceConfig(speaker="Guest", voice_config=types.VoiceConfig(
                    prebuilt_voice_config=types.PrebuiltVoiceConfig(voice_name="Puck")))
            ]
        )
    )
)
```

### 3.4 語音轉文字 (STT)

Gemini 模型原生支援音頻輸入，直接將音檔作為 content 傳入即可轉錄：

```python
import pathlib

audio_bytes = pathlib.Path("recording.mp3").read_bytes()

response = client.models.generate_content(
    model="gemini-2.5-flash",
    contents=[
        types.Content(parts=[
            types.Part(inline_data=types.Blob(mime_type="audio/mp3", data=audio_bytes)),
            types.Part(text="Transcribe this audio.")
        ])
    ]
)
print(response.text)
```

---

## 4. xAI (Grok)

| 項目 | 說明 |
|------|------|
| **API Base URL** | `https://api.x.ai/v1` |
| **認證方式** | `Authorization: Bearer xai-...` |
| **API 風格** | 相容 OpenAI Chat Completions 格式 |
| **Console** | `console.x.ai` |

### 4.1 LLM 模型

| 模型 ID | 類型 | 備註 |
|---------|------|------|
| `grok-4` | 旗艦 | 256K context，原生工具使用 |
| `grok-4.1-fast` | 快速版 | 4.1 更新版本 |
| `grok-3` | 前代旗艦 | GA |
| `grok-3-mini` | 輕量版 | 高性價比 |
| `grok-code-fast-1` | 程式碼 | 專為 agentic coding 優化 |

```python
from openai import OpenAI  # 使用 OpenAI SDK 相容介面

client = OpenAI(
    api_key="xai-...",
    base_url="https://api.x.ai/v1"
)

response = client.chat.completions.create(
    model="grok-4",
    messages=[
        {"role": "system", "content": "You are Grok, a helpful assistant."},
        {"role": "user", "content": "What's happening in the world today?"}
    ]
)
print(response.choices[0].message.content)
```

**Live Search (即時搜尋)**：
```python
response = client.chat.completions.create(
    model="grok-4",
    messages=[{"role": "user", "content": "Latest tech news today"}],
    tools=[
        {"type": "web_search"},    # 網頁搜尋
        {"type": "x_search"}       # X 平台搜尋
    ]
)
```

### 4.2 文生圖 / 影像

| 模型 / 端點 | 備註 |
|------------|------|
| Aurora (grok-2-image) | 圖片生成與編輯 |
| `grok-imagine-video` | 影片生成 (6-15秒) |

```python
# 圖片生成
response = client.images.generate(
    model="grok-2-image-1212",
    prompt="A photorealistic portrait in natural lighting"
)
```

### 4.3 語音 (Voice Agent API)

| 功能 | 端點 |
|------|------|
| Voice Agent (即時語音) | `wss://api.x.ai/v1/realtime` (WebSocket) |
| 獨立 TTS / STT | 即將推出 |

**Voice Agent 連線**：
```python
import websockets, json, os

XAI_API_KEY = os.getenv("XAI_API_KEY")

async with websockets.connect(
    uri="wss://api.x.ai/v1/realtime",
    additional_headers={"Authorization": f"Bearer {XAI_API_KEY}"}
) as ws:
    session_config = {
        "type": "session.update",
        "session": {
            "voice": "Ara",          # Ara / Eve / Leo 等
            "instructions": "You are a helpful assistant.",
            "tools": [{"type": "web_search"}]
        }
    }
    await ws.send(json.dumps(session_config))
```

> **可用語音**: Ara, Eve, Leo 等多種風格，支援 100+ 語言自動偵測。

---

## 5. ElevenLabs

| 項目 | 說明 |
|------|------|
| **API Base URL** | `https://api.elevenlabs.io/v1` |
| **認證方式** | `xi-api-key: ...` (Header) |
| **SDK** | `pip install elevenlabs` |
| **專長** | 業界頂尖 TTS / 語音克隆 / STT |

### 5.1 TTS 模型

| 模型 ID | 備註 |
|---------|------|
| `eleven_v3` | 最新，最具表現力，支援情感標籤 |
| `eleven_turbo_v2_5` | 高品質低延遲平衡 |
| `eleven_flash_v2_5` | 超低延遲 (75ms)，適合即時對話 |
| `eleven_multilingual_v2` | 多語言，32 種語言 |

```python
import requests

url = "https://api.elevenlabs.io/v1/text-to-speech/{voice_id}"

headers = {
    "xi-api-key": "YOUR_API_KEY",
    "Content-Type": "application/json"
}

data = {
    "text": "Hello, this is a test of ElevenLabs text to speech.",
    "model_id": "eleven_turbo_v2_5",
    "voice_settings": {
        "stability": 0.5,
        "similarity_boost": 0.75
    }
}

response = requests.post(url, json=data, headers=headers)
with open("output.mp3", "wb") as f:
    f.write(response.content)
```

**v3 Audio Tags (情感控制)**：
```
[whispers] Something's coming... [sighs] I can feel it.
[laughs] That's hilarious!
[excited] We just launched the new product!
```

**Text to Dialogue API (v3 多角色)**：
```python
data = {
    "model_id": "eleven_v3",
    "dialogue": [
        {"speaker": "narrator", "voice_id": "...", "text": "The room fell silent."},
        {"speaker": "alice", "voice_id": "...", "text": "[whispers] Did you hear that?"},
        {"speaker": "bob", "voice_id": "...", "text": "[nervous] I think we should leave."}
    ]
}
response = requests.post(
    "https://api.elevenlabs.io/v1/text-to-dialogue",
    json=data, headers=headers
)
```

### 5.2 語音轉文字 (Scribe)

| 模型 | 備註 |
|------|------|
| Scribe v2 | 支援 90+ 語言的高精度轉錄 |

```python
url = "https://api.elevenlabs.io/v1/speech-to-text"
headers = {"xi-api-key": "YOUR_API_KEY"}

with open("audio.mp3", "rb") as f:
    response = requests.post(
        url,
        headers=headers,
        files={"file": f},
        data={"model_id": "scribe_v2"}
    )
print(response.json()["text"])
```

### 5.3 其他功能

- **Voice Cloning**：上傳音頻樣本即可克隆聲音
- **Sound Effects**：文字生成音效
- **Music Generation**：文字生成音樂
- **Dubbing**：自動配音翻譯
- **Voice Agent Platform**：建構語音 AI 代理

---

## 6. WaveSpeed AI

| 項目 | 說明 |
|------|------|
| **API Base URL** | `https://api.wavespeed.ai/api/v3` |
| **認證方式** | `Authorization: Bearer YOUR_API_KEY` |
| **SDK** | Python / JavaScript / ComfyUI / N8N 整合 |
| **定位** | 多模態生成加速聚合平台 (700+ 模型) |
| **官方文件** | `wavespeed.ai/docs` |

> WaveSpeed AI 是一個**模型聚合平台**，透過單一 API 存取 700+ 來自不同供應商的圖像/影片/音頻生成模型 (FLUX, Kling, Veo, Seedance, WAN, Minimax 等)，主打極快推論速度（圖片 <2 秒，影片 <2 分鐘）。

### 6.1 API 工作流程 (非同步)

WaveSpeed 使用**提交任務 → 輪詢結果**的非同步模式：

**Step 1: 提交任務**
```python
import requests

response = requests.post(
    "https://api.wavespeed.ai/api/v3/wavespeed-ai/flux-dev",
    headers={
        "Authorization": "Bearer YOUR_API_KEY",
        "Content-Type": "application/json"
    },
    json={"prompt": "A cat wearing a space suit"}
)
data = response.json()
task_id = data["data"]["id"]
print(f"Task ID: {task_id}")
# 回應: {"code": 200, "data": {"id": "abc123", "status": "pending", "urls": {"get": "..."}}}
```

**Step 2: 輪詢取得結果**
```python
import time

while True:
    result = requests.get(
        f"https://api.wavespeed.ai/api/v3/predictions/{task_id}/result",
        headers={"Authorization": "Bearer YOUR_API_KEY"}
    ).json()

    if result["data"]["status"] == "completed":
        print("Done!", result["data"]["outputs"])
        break
    elif result["data"]["status"] == "failed":
        print("Failed:", result["data"]["error"])
        break
    time.sleep(1)
```

### 6.2 常用模型端點

模型端點格式：`https://api.wavespeed.ai/api/v3/{provider}/{model-name}`

| 類型 | 模型端點範例 | 說明 |
|------|-------------|------|
| **圖片生成** | `wavespeed-ai/flux-dev` | FLUX Dev 快速圖片 |
| **圖片生成** | `bytedance/seedream-4.5` | ByteDance Seedream |
| **圖片生成** | `google/gemini-3-pro-image` | Nano Banana Pro |
| **圖片生成** | `openai/gpt-image-1.5` | OpenAI GPT Image 1.5 |
| **影片生成** | `kling/omni-video-o3` | Kling Omni V3 |
| **影片生成** | `alibaba/wan-2.6-i2v-pro` | WAN 2.6 圖轉影片 |
| **影片生成** | `bytedance/seedance` | ByteDance Seedance |
| **影片生成** | `google/veo-3.1` | Google Veo 3.1 |
| **語音生成** | `minimax/speech-02` | Minimax Speech 02 |
| **音樂生成** | `ace-step/ace-step-1.5` | ACE-Step 音樂生成 |
| **影片升級** | `wavespeed-ai/video-upscaler` | 4K 影片升級 |

**Python SDK 用法**：
```python
import wavespeed

output = wavespeed.run(
    "openai/gpt-image-1.5",
    {"prompt": "A serene Japanese garden at sunset"}
)
print(output["outputs"][0])  # 圖片 URL
```

### 6.3 重點特性

- **Webhook 支援**：可設定 callback URL 接收完成通知
- **LoRA 訓練**：支援自訂 LoRA 模型訓練與使用
- **MCP 整合**：支援 AI Agent 透過 MCP 協議呼叫生成功能
- **ComfyUI 整合**：提供 ComfyUI 自訂節點
- **輸出暫存**：生成結果保存 7 天

---

## 7. fal.ai

| 項目 | 說明 |
|------|------|
| **API Base URL** | `https://queue.fal.run/{model-id}` |
| **認證方式** | `Authorization: Key YOUR_FAL_KEY` |
| **SDK** | `pip install fal-client` / `npm install @fal-ai/client` |
| **定位** | 生成式媒體 API 平台 (600+ 模型)，速度極快 |
| **官方文件** | `docs.fal.ai` |

> fal.ai 提供 600+ 生產級生成媒體模型，專注於**極速推論**，號稱全球最快的 FLUX 推論引擎。支援圖片、影片、語音、3D 等多模態生成。

### 7.1 基本使用 (Python)

```python
import fal_client
import os

os.environ["FAL_KEY"] = "YOUR_FAL_KEY"

# 同步模式 (subscribe 會自動等待完成)
result = fal_client.subscribe(
    "fal-ai/flux/dev",
    arguments={
        "prompt": "Photo of a rhino in a suit sitting at a bar, award winning photography",
        "image_size": "landscape_4_3",
        "num_images": 1,
    }
)
print(result["images"][0]["url"])
```

### 7.2 JavaScript 使用

```javascript
import { fal } from "@fal-ai/client";

const result = await fal.subscribe("fal-ai/flux/dev", {
    input: {
        prompt: "A beautiful sunset over mountains",
        image_size: "landscape_16_9"
    },
    logs: true,
    onQueueUpdate: (update) => {
        if (update.status === "IN_PROGRESS") {
            update.logs.map((log) => log.message).forEach(console.log);
        }
    },
});
console.log(result.data.images[0].url);
```

### 7.3 Queue API (非同步模式)

```python
# 提交任務
request = fal_client.submit(
    "fal-ai/flux/dev",
    arguments={"prompt": "A futuristic cityscape"}
)
request_id = request.request_id

# 查詢狀態
status = fal_client.status("fal-ai/flux/dev", request_id, with_logs=True)

# 取得結果
result = fal_client.result("fal-ai/flux/dev", request_id)
```

### 7.4 常用模型 ID

| 類型 | 模型 ID | 說明 |
|------|---------|------|
| **圖片 (FLUX)** | `fal-ai/flux/dev` | FLUX.1 Dev (通用) |
| **圖片 (FLUX)** | `fal-ai/flux/schnell` | FLUX.1 Schnell (超快 1-4 步) |
| **圖片 (FLUX)** | `fal-ai/flux-pro/v1.1` | FLUX.1 Pro 1.1 (最高品質) |
| **圖片 (FLUX 2)** | `fal-ai/flux-2-max` | FLUX.2 MAX (最新旗艦) |
| **圖片 (FLUX 2)** | `fal-ai/flux-2/flash` | FLUX.2 Flash (快速) |
| **圖片編輯** | `fal-ai/flux/dev/image-to-image` | FLUX 圖轉圖 |
| **圖片+LoRA** | `fal-ai/flux-lora` | FLUX + 自訂 LoRA |
| **圖片編輯** | `fal-ai/flux-kontext` | FLUX Kontext (上下文編輯) |
| **影片** | `fal-ai/kling-video/v2.1/standard` | Kling v2.1 |
| **影片** | `fal-ai/minimax/video-01` | MiniMax Hailuo |
| **影片** | `xai/grok-imagine-video` | Grok Imagine 圖轉影片 |
| **語音 STT** | `fal-ai/whisper` | Whisper 語音轉文字 |

### 7.5 重點特性

- **Webhook 支援**：提交任務時附帶 `webhookUrl` 參數
- **檔案上傳**：`fal_client.upload(file)` 或 `fal.storage.upload(file)`
- **串流**：支援 `fal.stream()` 即時串流生成
- **Cursor MCP 整合**：可直接在 Cursor IDE 中使用 fal 模型
- **LoRA 訓練**：支援線上 LoRA 微調

---

## 8. Together AI

| 項目 | 說明 |
|------|------|
| **API Base URL** | `https://api.together.xyz/v1` |
| **認證方式** | `Authorization: Bearer YOUR_TOGETHER_API_KEY` |
| **SDK** | `pip install together` / `npm install together-ai` |
| **定位** | 開源模型推論平台 (200+ 模型)，OpenAI API 完全相容 |
| **官方文件** | `docs.together.ai` |

> Together AI 是**開源模型推論平台**，提供 Llama、Qwen、DeepSeek、Mistral 等 200+ 開源 LLM，以及圖片/影片生成模型。最大亮點是 **OpenAI SDK 完全相容**，可無縫切換。

### 8.1 基本使用 (OpenAI SDK 相容)

```python
from openai import OpenAI
import os

client = OpenAI(
    api_key=os.environ.get("TOGETHER_API_KEY"),
    base_url="https://api.together.xyz/v1",
)

response = client.chat.completions.create(
    model="meta-llama/Llama-4-Maverick-17B-128E-Instruct-FP8",
    messages=[
        {"role": "system", "content": "You are a helpful assistant."},
        {"role": "user", "content": "Tell me the top 3 things to do in Tokyo"},
    ],
)
print(response.choices[0].message.content)
```

### 8.2 Together 原生 SDK

```python
from together import Together

client = Together()  # 自動讀取 TOGETHER_API_KEY 環境變數

response = client.chat.completions.create(
    model="meta-llama/Llama-Vision-Free",
    messages=[{"role": "user", "content": "What are some fun things to do?"}],
)
print(response.choices[0].message.content)
```

### 8.3 圖片生成

```python
response = client.images.generate(
    model="black-forest-labs/FLUX.1-schnell-Free",
    prompt="A beautiful mountain landscape at golden hour",
    n=1,
)
print(response.data[0].url)
```

### 8.4 Vision (視覺理解)

```python
response = client.chat.completions.create(
    model="meta-llama/Llama-4-Maverick-17B-128E-Instruct-FP8",
    messages=[{
        "role": "user",
        "content": [
            {"type": "text", "text": "What's in this image?"},
            {"type": "image_url", "image_url": {"url": "https://example.com/photo.jpg"}}
        ]
    }]
)
```

### 8.5 常用模型 ID

| 類型 | 模型 ID | 說明 |
|------|---------|------|
| **LLM (旗艦)** | `meta-llama/Llama-4-Maverick-17B-128E-Instruct-FP8` | Llama 4 Maverick |
| **LLM (推理)** | `deepseek-ai/DeepSeek-R1` | DeepSeek R1 推理模型 |
| **LLM (輕量)** | `Qwen/Qwen3-Next-80B-A3B-Instruct` | Qwen3 Next MoE |
| **LLM (開放)** | `openai/gpt-oss-20b` | OpenAI GPT-OSS 20B |
| **LLM (免費)** | `meta-llama/Llama-Vision-Free` | Llama Vision (免費) |
| **圖片** | `black-forest-labs/FLUX.1-schnell-Free` | FLUX Schnell (免費) |
| **圖片** | `black-forest-labs/FLUX.1.1-pro` | FLUX 1.1 Pro |
| **影片** | 透過 Runware 合作提供 Veo, Sora 2, Seedance 等 | 新增功能 |

### 8.6 重點特性

- **OpenAI SDK 完全相容**：只需改 `base_url` 和 `api_key` 即可從 OpenAI 遷移
- **開源模型為主**：Llama 4, DeepSeek R1, Qwen3, Mistral 等
- **Fine-tuning**：支援 LoRA 微調開源模型
- **Async 支援**：`AsyncTogether` 客戶端支援並行請求
- **Tool Calling**：支援 function calling / tool use
- **Embedding**：支援文字嵌入向量模型

---

## 9. 跨平台功能對照表

| 功能 | OpenAI | Anthropic | Google Gemini | xAI Grok | ElevenLabs | WaveSpeed | fal.ai | Together AI |
|------|--------|-----------|--------------|----------|------------|-----------|--------|-------------|
| **LLM** | ✅ GPT-5.2 | ✅ Opus 4.6 | ✅ Gemini 3.1 Pro | ✅ Grok 4 | ❌ | ⚠️ LLM 服務 | ❌ | ✅ Llama 4/DeepSeek R1 |
| **文生圖** | ✅ GPT-Image-1.5 | ❌ | ✅ Imagen 4 | ✅ Aurora | ❌ | ✅ 700+ 模型 | ✅ FLUX 2/600+ | ✅ FLUX |
| **TTS** | ✅ gpt-4o-mini-tts | ❌ | ✅ Gemini TTS | ⚠️ Voice Agent | ✅ Eleven v3 | ✅ 轉接模型 | ❌ | ❌ |
| **STT** | ✅ transcribe | ❌ | ✅ 原生音頻 | ⚠️ Voice Agent | ✅ Scribe v2 | ❌ | ✅ Whisper | ❌ |
| **影片生成** | ✅ Sora 2 | ❌ | ✅ Veo 3.1 | ✅ Grok Imagine | ❌ | ✅ 多模型 | ✅ Kling/Minimax | ✅ 多模型 |
| **OpenAI SDK 相容** | ✅ 原生 | ❌ | ❌ | ✅ | ❌ | ❌ | ❌ | ✅ 完全相容 |
| **開源模型** | ❌ | ❌ | ❌ | ❌ | ❌ | ✅ 聚合 | ✅ 聚合 | ✅ 專精 |
| **免費方案** | ❌ | ❌ | ✅ 有限 | ❌ | ✅ 有限 | ❌ | ✅ 有限 | ✅ 有限 |

### 快速選型建議

| 需求場景 | 推薦方案 |
|---------|---------|
| 最強推理能力 | Claude Opus 4.6 或 Gemini 3.1 Pro |
| 最佳程式碼生成 | Claude Sonnet 4.5 或 GPT-5.2 |
| 最高品質圖片 | GPT-Image-1.5 或 Imagen 4 Ultra |
| 最快圖片生成 | fal.ai FLUX Schnell 或 WaveSpeed FLUX |
| 最多圖片/影片模型選擇 | WaveSpeed AI (700+) 或 fal.ai (600+) |
| 最佳 TTS 品質 | ElevenLabs Eleven v3 |
| 最低延遲 TTS | ElevenLabs Flash v2.5 (75ms) |
| 最佳 STT | OpenAI gpt-4o-mini-transcribe |
| 即時語音代理 | xAI Voice Agent (最低延遲 <1s) |
| 低成本高速 LLM | Claude Haiku 4.5 或 Gemini 2.5 Flash |
| 開源 LLM 推論 | Together AI (Llama 4, DeepSeek R1) |
| 全功能一站式 | OpenAI (LLM + 圖 + 音頻 + 影片) |
| 影片生成聚合 | WaveSpeed (Kling/Veo/Seedance/WAN) |

---

## 環境變數設定建議

```bash
# .env 檔案範例
OPENAI_API_KEY=sk-...
ANTHROPIC_API_KEY=sk-ant-...
GOOGLE_API_KEY=AIza...
XAI_API_KEY=xai-...
ELEVENLABS_API_KEY=...
WAVESPEED_API_KEY=...
FAL_KEY=...
TOGETHER_API_KEY=...
```

```python
# Claude Code 中統一載入
import os
from dotenv import load_dotenv
load_dotenv()

OPENAI_KEY = os.getenv("OPENAI_API_KEY")
ANTHROPIC_KEY = os.getenv("ANTHROPIC_API_KEY")
GOOGLE_KEY = os.getenv("GOOGLE_API_KEY")
XAI_KEY = os.getenv("XAI_API_KEY")
ELEVENLABS_KEY = os.getenv("ELEVENLABS_API_KEY")
WAVESPEED_KEY = os.getenv("WAVESPEED_API_KEY")
FAL_KEY = os.getenv("FAL_KEY")
TOGETHER_KEY = os.getenv("TOGETHER_API_KEY")
```

---

*此文件僅供參考，各平台 API 可能隨時更新，建議以官方文件為準。*
