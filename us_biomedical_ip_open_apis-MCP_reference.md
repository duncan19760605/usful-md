# 美國生醫/智財開放 API 參考手冊

> 本文件彙整 FDA、USPTO、NIH、PubMed、PubChem 等機構提供的開放 API，供 Claude Code 或自動化腳本直接查詢使用。
>
> 最後更新：2026-02-20

---

## 目錄

1. [總覽比較表](#1-總覽比較表)
2. [openFDA API（藥物 / 醫材 / 食品 / 菸草）](#2-openfda-api)
3. [USPTO API（專利 / 商標）](#3-uspto-api)
4. [NIH API（研究經費 / 臨床試驗）](#4-nih-api)
5. [NCBI E-utilities / PubMed API（文獻）](#5-ncbi-e-utilities--pubmed-api)
6. [PubChem API（化合物 / 生物活性）](#6-pubchem-api)
7. [ClinicalTrials.gov API（臨床試驗）](#7-clinicaltrialsgov-api)
8. [實用程式碼範例](#8-實用程式碼範例)
9. [速率限制與最佳實務](#9-速率限制與最佳實務)
10. [**MCP Server 整合指南**](#10-mcp-server-整合指南)（Claude Code / Claude Desktop）

---

## 1. 總覽比較表

| 機構/服務 | 領域 | 費用 | 需要 API Key | 註冊方式 | 回傳格式 | Base URL |
|---|---|---|---|---|---|---|
| **openFDA** | 藥物、醫材、食品、菸草 | 免費 | 建議（可提升額度） | <https://open.fda.gov/apis/authentication/> | JSON | `https://api.fda.gov/` |
| **USPTO PatentsView** | 專利分析 | 免費 | 不需要 | — | JSON | `https://search.patentsview.org/api/v1/` |
| **USPTO Open Data Portal (ODP)** | 專利申請案/審查資料 | 免費 | 需要 | <https://developer.uspto.gov> | JSON/XML | `https://beta-api.uspto.gov/api/v1/` |
| **USPTO TSDR** | 商標狀態 | 免費 | 不需要 | — | XML | `https://tsdrapi.uspto.gov/` |
| **NIH RePORTER** | 研究經費/補助 | 免費 | 不需要 | — | JSON | `https://api.reporter.nih.gov/v2/` |
| **NCBI E-utilities (PubMed)** | 生醫文獻 | 免費 | 建議（可提升額度） | NCBI 帳號 → API Key | XML/JSON | `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/` |
| **PubChem PUG-REST** | 化合物、物質、生物活性 | 免費 | 不需要 | — | JSON/XML/CSV/SDF | `https://pubchem.ncbi.nlm.nih.gov/rest/pug/` |
| **ClinicalTrials.gov v2** | 臨床試驗 | 免費 | 不需要 | — | JSON/CSV | `https://clinicaltrials.gov/api/v2/` |

---

## 2. openFDA API

### 2.1 簡介

openFDA 是 FDA 提供的公開資料 API，以 Elasticsearch 為底層，涵蓋藥品、醫療器材、食品、菸草等產品的不良事件、召回、標示及核准資訊。

### 2.2 所有端點一覽

#### 藥物（Drug）
| 端點 | 用途 | URL |
|---|---|---|
| `/drug/event.json` | 藥物不良事件報告 (FAERS) | 查詢副作用、用藥錯誤、治療失敗 |
| `/drug/label.json` | 藥品標示 (SPL) | 仿單內容、適應症、禁忌、用法用量 |
| `/drug/ndc.json` | NDC 目錄 | 國家藥品編碼、製造商、包裝資訊 |
| `/drug/drugsfda.json` | Drugs@FDA | 已核准藥品、核准歷程、治療等效碼 |
| `/drug/enforcement.json` | 藥品召回 | 召回分類、原因、範圍 |

#### 醫療器材（Device）
| 端點 | 用途 |
|---|---|
| `/device/event.json` | 醫材不良事件 (MDR) |
| `/device/classification.json` | 醫材分類（Class I/II/III） |
| `/device/510k.json` | 510(k) 上市前通知 |
| `/device/pma.json` | PMA 上市前核准 |
| `/device/recall.json` | 醫材召回 |
| `/device/registrationlisting.json` | 場所登記與器材列表 |
| `/device/enforcement.json` | 醫材強制召回 |
| `/device/udi.json` | UDI 唯一器材識別 |
| `/device/covid19serology.json` | COVID-19 血清學檢測 |

#### 食品（Food）
| 端點 | 用途 |
|---|---|
| `/food/event.json` | 食品不良事件 (CAERS) |
| `/food/enforcement.json` | 食品召回 |

#### 菸草（Tobacco）
| 端點 | 用途 |
|---|---|
| `/tobacco/problem.json` | 菸草產品問題報告 |

#### 其他（Other）
| 端點 | 用途 |
|---|---|
| `/other/substance.json` | 物質資料 (UNII) |
| `/other/nsde.json` | 藥品標示新增/修改事件 |

### 2.3 認證與額度

- **無 API Key**：每天 1,000 次請求，每分鐘 40 次
- **有 API Key**：每天 120,000 次請求，每分鐘 240 次
- **取得方式**：至 <https://open.fda.gov/apis/authentication/> 免費註冊即得
- 單次查詢最大回傳 1,000 筆（`limit=1000`）

### 2.4 查詢語法

```
https://api.fda.gov/{category}/{endpoint}.json?search={field}:{term}&limit={n}&api_key={key}
```

常用參數：
- `search=` — 搜尋條件，支援 `AND`、`OR`、`NOT`、萬用字元
- `count=` — 按欄位計數（取代 search 結果，改為回傳彙總統計）
- `limit=` — 回傳筆數上限（最大 1000）
- `skip=` — 分頁偏移量（最大 26,000）

### 2.5 範例

```bash
# 查詢 aspirin 相關不良事件（取 5 筆）
curl "https://api.fda.gov/drug/event.json?search=patient.drug.openfda.generic_name:aspirin&limit=5"

# 統計某醫材的不良事件反應類型
curl "https://api.fda.gov/device/event.json?count=event_type.exact"

# 搜尋特定 510(k) 核准
curl "https://api.fda.gov/device/510k.json?search=applicant:medtronic&limit=10"
```

```python
import requests

# Python 範例：查詢藥品標示
url = "https://api.fda.gov/drug/label.json"
params = {
    "search": "openfda.brand_name:tylenol",
    "limit": 3,
    "api_key": "YOUR_API_KEY"  # 可選
}
resp = requests.get(url, params=params)
data = resp.json()
for r in data["results"]:
    print(r["openfda"].get("generic_name", ["N/A"]))
```

---

## 3. USPTO API

### 3.1 專利相關 API

#### 3.1.1 PatentsView API（專利分析）

- **用途**：搜尋自 1976 年起的已授權美國專利，支援發明人、受讓人、CPC 分類、引用等多維度查詢
- **費用**：完全免費
- **API Key**：不需要
- **Base URL**：`https://search.patentsview.org/api/v1/`
- **端點**：`patent`, `inventor`, `assignee`, `cpc_subclass`, `uspc`, `location`
- **文件**：<https://patentsview.org/apis/purpose>

```bash
# 查詢含 "3D printing" 的專利
curl "https://search.patentsview.org/api/v1/patent/?q={\"_text_any\":{\"patent_abstract\":\"3D printing\"}}&f=[\"patent_id\",\"patent_title\",\"patent_date\"]&o={\"per_page\":5}"
```

#### 3.1.2 USPTO Open Data Portal (ODP)（專利審查資料）

- **用途**：取得專利申請書目資料、審查歷程、文件下載、期限調整、轉讓資訊
- **費用**：免費
- **API Key**：需要（免費申請）
- **申請方式**：至 USPTO Developer Hub <https://developer.uspto.gov> 註冊帳號取得 API Key
- **涵蓋資料**：1981 年至今已公開的專利申請案（部分追溯至 1935 年）
- **更新頻率**：每日

主要端點：
| 端點 | 用途 |
|---|---|
| `/search` | 搜尋專利申請 |
| `/meta-data/{appId}` | 取得申請案書目資料 |
| `/documents/{appId}` | 下載審查歷程文件（IFW） |
| `/assignments/{appId}` | 轉讓/讓與資訊 |
| `/term-adjustments/{appId}` | 專利期限調整 (PTA) |
| `/transactions/{appId}` | 交易紀錄 |
| `/foreign-priority/{appId}` | 優先權資料 |
| `/continuity/{appId}` | 延續案資料 |
| `/bulk-data` | 大量資料下載 |

```python
import requests

headers = {"X-API-Key": "YOUR_USPTO_API_KEY"}

# 搜尋專利申請
resp = requests.get(
    "https://beta-api.uspto.gov/api/v1/search",
    params={"q": "inventionTitle:artificial intelligence", "rows": 5},
    headers=headers
)
print(resp.json())
```

#### 3.1.3 PTAB API（專利審判委員會）

- **用途**：IPR（多方複審）、PGR（授權後複審）、申訴案的裁決文件及資料
- **費用**：免費
- **文件**：包含在 ODP 中

### 3.2 商標相關 API

#### 3.2.1 TSDR API（商標狀態查詢）

- **用途**：查詢商標申請案的狀態、所有人、商品/服務描述
- **費用**：免費
- **API Key**：不需要
- **URL 格式**：`https://tsdrapi.uspto.gov/ts/cd/casestatus/{serialNumber}/info`

```bash
# 以序號查詢商標狀態
curl "https://tsdrapi.uspto.gov/ts/cd/casestatus/87654321/info"
```

#### 3.2.2 PatentsView 商標

> 注意：USPTO 官方的商標全文搜尋 API 資源較有限。常見替代方案：
> - **TESS (Trademark Electronic Search System)**：網頁版全文搜尋（無正式 REST API）
> - **RapidAPI 上的非官方 USPTO Trademark API**：第三方服務，有免費/付費方案
> - **Marker API**：付費商標搜尋 API

### 3.3 大量資料下載

USPTO 提供大量專利與商標資料的 bulk data 下載：
- **專利**：<https://developer.uspto.gov/data>（XML/JSON 格式，自 1976 年起）
- **商標**：<https://developer.uspto.gov/data>（每日/每週更新）
- **格式**：XML, JSON, TSV
- **授權**：公共領域（Public Domain），建議註明出處

---

## 4. NIH API

### 4.1 NIH RePORTER API（研究經費）

- **用途**：查詢 NIH 及其他聯邦機構的研究補助計畫、金額、PI 資訊、出版物、專利
- **費用**：完全免費
- **API Key**：不需要
- **Base URL**：`https://api.reporter.nih.gov/v2/`
- **文件**：<https://api.reporter.nih.gov/>

主要端點：
| 端點 | 方法 | 用途 |
|---|---|---|
| `/projects/search` | POST | 搜尋研究計畫 |
| `/publications/search` | POST | 搜尋關聯出版物 |

```python
import requests

# 搜尋 NIH 補助計畫
url = "https://api.reporter.nih.gov/v2/projects/search"
payload = {
    "criteria": {
        "project_terms": "CRISPR gene therapy",
        "fiscal_years": [2024, 2025],
        "include_active_projects": True
    },
    "offset": 0,
    "limit": 10
}
resp = requests.post(url, json=payload)
data = resp.json()
for p in data.get("results", []):
    print(f"{p['project_num']}: {p['project_title']} — ${p.get('award_amount', 'N/A')}")
```

使用建議：
- 每秒不超過 1 次請求
- 大量作業建議安排於週末或工作日 21:00–05:00 EST

---

## 5. NCBI E-utilities / PubMed API

### 5.1 簡介

NCBI E-utilities 是存取 Entrez 系統中所有 38 個資料庫的公共 API，涵蓋 PubMed、PMC、Gene、Protein、Nucleotide 等。

### 5.2 核心工具

| 工具 | 用途 | URL 模式 |
|---|---|---|
| **EInfo** | 列出資料庫資訊/索引欄位 | `einfo.fcgi?db={db}` |
| **ESearch** | 搜尋，回傳 UID 列表 | `esearch.fcgi?db={db}&term={query}` |
| **ESummary** | 取得文獻摘要 | `esummary.fcgi?db={db}&id={uids}` |
| **EFetch** | 取得完整紀錄 | `efetch.fcgi?db={db}&id={uids}&rettype={type}` |
| **ELink** | 跨資料庫連結 | `elink.fcgi?dbfrom={db}&db={db2}&id={uids}` |
| **EPost** | 上傳 UID 列表至 History Server | `epost.fcgi?db={db}&id={uids}` |
| **ESpell** | 拼字建議 | `espell.fcgi?db={db}&term={query}` |
| **EGQuery** | 全資料庫搜尋計數 | `egquery.fcgi?term={query}` |
| **ECitMatch** | 引文匹配取得 PMID | `ecitmatch.cgi?...` |

### 5.3 認證與額度

- **無 API Key**：每秒 3 次請求，每日無硬性上限但可能被限速
- **有 API Key**：每秒 10 次請求
- **取得方式**：
  1. 至 <https://www.ncbi.nlm.nih.gov/account/> 建立 NCBI 帳號（需透過第三方登入如 Google/ORCID）
  2. 登入後進入 Settings → API Key Management → Create an API Key
  3. 在請求中加入 `&api_key=YOUR_KEY`

### 5.4 常用範例

```bash
# 搜尋 PubMed 中關於 "CRISPR delivery" 的文獻
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=CRISPR+delivery&retmax=5&retmode=json&api_key=YOUR_KEY"

# 取得特定 PMID 的摘要（XML 格式）
curl "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=33456789&rettype=abstract&retmode=xml"
```

```python
from Bio import Entrez

Entrez.email = "your@email.com"
Entrez.api_key = "YOUR_NCBI_API_KEY"

# 搜尋 PubMed
handle = Entrez.esearch(db="pubmed", term="mRNA vaccine delivery", retmax=10)
record = Entrez.read(handle)
ids = record["IdList"]

# 取得摘要
handle = Entrez.efetch(db="pubmed", id=ids, rettype="medline", retmode="text")
print(handle.read())
```

### 5.5 Entrez Direct (EDirect) — 命令列工具

```bash
# 安裝
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

# 搜尋並取得結果
esearch -db pubmed -query "CAR-T cell therapy" | efetch -format abstract | head -50
```

---

## 6. PubChem API

### 6.1 簡介

PubChem 提供多種程式化存取方式，其中 PUG-REST 是最常用的 RESTful API，可查詢化合物、物質及生物活性資料。

### 6.2 PUG-REST

- **費用**：完全免費
- **API Key**：不需要
- **Base URL**：`https://pubchem.ncbi.nlm.nih.gov/rest/pug/`

URL 結構：
```
https://pubchem.ncbi.nlm.nih.gov/rest/pug/{domain}/{namespace}/{identifiers}/{operation}/{output}
```

| 組件 | 選項 |
|---|---|
| **domain** | `compound`, `substance`, `assay`, `gene`, `protein`, `pathway`, `cell`, `taxonomy` |
| **namespace** | `cid`, `sid`, `aid`, `name`, `smiles`, `inchikey`, `formula`, `substructure`, `similarity` |
| **operation** | `record`, `property`, `synonyms`, `sids`, `cids`, `aids`, `description`, `conformers`, `xrefs` |
| **output** | `JSON`, `XML`, `CSV`, `TXT`, `SDF`, `PNG` |

### 6.3 PUG-View

用於取得 PubChem 紀錄頁面中的註釋資訊（如安全性、毒理學等第三方資料）。

```
https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{CID}/JSON
```

### 6.4 範例

```bash
# 以名稱查詢化合物屬性
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/property/MolecularFormula,MolecularWeight,IUPACName/JSON"

# 以 SMILES 搜尋相似化合物
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/CC(=O)OC1=CC=CC=C1C(O)=O/JSON?Threshold=90&MaxRecords=5"

# 取得化合物 2D 結構圖
curl -o aspirin.png "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/aspirin/PNG"

# 查詢生物活性資料
curl "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/assaysummary/JSON"
```

```python
import pubchempy as pcp

# 以名稱搜尋
compounds = pcp.get_compounds("ibuprofen", "name")
for c in compounds:
    print(f"CID: {c.cid}, Formula: {c.molecular_formula}, MW: {c.molecular_weight}")

# 以 SMILES 搜尋
compounds = pcp.get_compounds("CC(=O)OC1=CC=CC=C1C(O)=O", "smiles")
```

### 6.5 速率限制

- 每秒不超過 5 次請求
- 不適合大量批次下載（請使用 PubChem FTP）
- 加入 `time.sleep(0.2)` 控制速率

---

## 7. ClinicalTrials.gov API

### 7.1 簡介

ClinicalTrials.gov API v2 是現代化的 REST API（OpenAPI 3.0），提供全球臨床試驗資料的程式化存取。

### 7.2 基本資訊

- **費用**：完全免費
- **API Key**：不需要
- **Base URL**：`https://clinicaltrials.gov/api/v2/`
- **文件**：<https://clinicaltrials.gov/data-api/api>

### 7.3 端點

| 端點 | 用途 |
|---|---|
| `/studies` | 搜尋臨床試驗 |
| `/studies/{nctId}` | 取得特定試驗的詳細資訊 |
| `/stats/size` | 取得總研究數量 |
| `/stats/fieldValues` | 取得特定欄位的可用值 |
| `/version` | API 版本資訊 |

### 7.4 範例

```bash
# 搜尋正在招募的肺癌免疫治療臨床試驗
curl "https://clinicaltrials.gov/api/v2/studies?query.cond=lung+cancer&query.intr=immunotherapy&filter.overallStatus=RECRUITING&pageSize=5&format=json"

# 取得特定試驗
curl "https://clinicaltrials.gov/api/v2/studies/NCT04267848?format=json"
```

```python
import requests

# 搜尋臨床試驗
resp = requests.get("https://clinicaltrials.gov/api/v2/studies", params={
    "query.cond": "diabetes",
    "query.intr": "GLP-1",
    "filter.overallStatus": "RECRUITING",
    "pageSize": 10,
    "format": "json",
    "countTotal": "true"
})
data = resp.json()
print(f"共 {data.get('totalCount', 'N/A')} 項試驗")
for study in data.get("studies", []):
    proto = study["protocolSection"]["identificationModule"]
    print(f"  {proto['nctId']}: {proto['briefTitle']}")
```

---

## 8. 實用程式碼範例

### 8.1 跨 API 整合查詢：從藥名到完整情報

```python
"""
從藥物名稱出發，整合多個 API 取得完整藥物情報
"""
import requests
import time

DRUG_NAME = "pembrolizumab"

# === Step 1: openFDA — 取得核准資訊與標示 ===
fda_resp = requests.get("https://api.fda.gov/drug/drugsfda.json", params={
    "search": f"openfda.generic_name:{DRUG_NAME}",
    "limit": 1
})
fda_data = fda_resp.json().get("results", [{}])[0]
print(f"[FDA] Brand: {fda_data.get('openfda', {}).get('brand_name', ['N/A'])}")

time.sleep(0.5)

# === Step 2: openFDA — 不良事件統計 ===
ae_resp = requests.get("https://api.fda.gov/drug/event.json", params={
    "search": f"patient.drug.openfda.generic_name:{DRUG_NAME}",
    "count": "patient.reaction.reactionmeddrapt.exact",
    "limit": 5
})
print("[FDA AE] Top reactions:")
for item in ae_resp.json().get("results", []):
    print(f"  {item['term']}: {item['count']}")

time.sleep(0.5)

# === Step 3: PubChem — 化合物資訊 ===
pc_resp = requests.get(
    f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{DRUG_NAME}/property/MolecularFormula,MolecularWeight/JSON"
)
if pc_resp.status_code == 200:
    props = pc_resp.json()["PropertyTable"]["Properties"][0]
    print(f"[PubChem] Formula: {props['MolecularFormula']}, MW: {props['MolecularWeight']}")

time.sleep(0.5)

# === Step 4: PubMed — 最新文獻 ===
pm_resp = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", params={
    "db": "pubmed",
    "term": f"{DRUG_NAME}[Title] AND 2025[PDAT]",
    "retmax": 5,
    "retmode": "json"
})
pmids = pm_resp.json()["esearchresult"]["idlist"]
print(f"[PubMed] 2025 年相關文獻 PMIDs: {pmids}")

time.sleep(0.5)

# === Step 5: ClinicalTrials.gov — 進行中的試驗 ===
ct_resp = requests.get("https://clinicaltrials.gov/api/v2/studies", params={
    "query.intr": DRUG_NAME,
    "filter.overallStatus": "RECRUITING",
    "pageSize": 5,
    "format": "json",
    "countTotal": "true"
})
ct_data = ct_resp.json()
print(f"[ClinicalTrials] 招募中試驗數: {ct_data.get('totalCount', 'N/A')}")

time.sleep(0.5)

# === Step 6: NIH RePORTER — 相關補助計畫 ===
nih_resp = requests.post("https://api.reporter.nih.gov/v2/projects/search", json={
    "criteria": {
        "project_terms": DRUG_NAME,
        "include_active_projects": True
    },
    "offset": 0,
    "limit": 3
})
for p in nih_resp.json().get("results", []):
    print(f"[NIH] {p['project_num']}: {p['project_title'][:80]}...")
```

### 8.2 Node.js 範例

```javascript
// 使用 node-fetch（v3+ 需 ESM）
const fetch = (...args) => import('node-fetch').then(({default: f}) => f(...args));

// openFDA 查詢
async function searchFDA(drugName) {
  const url = `https://api.fda.gov/drug/label.json?search=openfda.generic_name:${drugName}&limit=3`;
  const resp = await fetch(url);
  const data = await resp.json();
  return data.results;
}

// PubChem 查詢
async function searchPubChem(name) {
  const url = `https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/${name}/property/MolecularFormula,MolecularWeight,IUPACName/JSON`;
  const resp = await fetch(url);
  return await resp.json();
}

// PubMed 查詢
async function searchPubMed(query) {
  const url = `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&term=${encodeURIComponent(query)}&retmax=5&retmode=json`;
  const resp = await fetch(url);
  return await resp.json();
}
```

---

## 9. 速率限制與最佳實務

### 速率限制總整理

| API | 無 Key 限制 | 有 Key 限制 |
|---|---|---|
| openFDA | 40/min, 1,000/day | 240/min, 120,000/day |
| NCBI E-utilities | 3/sec | 10/sec |
| PubChem PUG-REST | 5/sec（建議） | 同左（無 Key 機制） |
| USPTO ODP | 有限制（未公開具體數字） | 較寬鬆 |
| NIH RePORTER | 1/sec（建議） | 同左（無 Key 機制） |
| ClinicalTrials.gov v2 | 無官方硬性限制 | N/A |

### 最佳實務

1. **務必加入延遲**：在迴圈中使用 `time.sleep()` 或等效機制，避免被封鎖
2. **設定 User-Agent / email**：NCBI 要求在請求中附帶 `tool` 和 `email` 參數
3. **使用 Entrez History**：大量 PubMed 查詢時，用 `usehistory=y` + `WebEnv` + `query_key` 分批下載
4. **偏好 JSON 格式**：多數 API 支援 JSON，解析更方便
5. **錯誤處理**：實作指數退避重試機制（exponential backoff）
6. **大量資料**：優先考慮 bulk data 下載而非逐筆 API 呼叫
   - openFDA: <https://open.fda.gov/data/downloads/>
   - PubChem: <https://ftp.ncbi.nlm.nih.gov/pubchem/>
   - USPTO: <https://developer.uspto.gov/data>
   - NIH ExPORTER: <https://exporter.nih.gov/>
7. **快取結果**：避免重複查詢相同資料

### Claude Code 使用建議

在 Claude Code 中使用這些 API 時：
- 可直接用 `curl` 或 Python `requests` 呼叫，大部分 API 不需要特別安裝套件
- 若需大量操作 PubMed，建議安裝 `biopython`（`pip install biopython`）
- 若需操作 PubChem，建議安裝 `pubchempy`（`pip install pubchempy`）
- 所有 API 均為 HTTPS，請確保網路存取已啟用
- API Key 可透過環境變數設定：
  ```bash
  export OPENFDA_API_KEY="your_key"
  export NCBI_API_KEY="your_key"
  export USPTO_API_KEY="your_key"
  ```

---

## 10. MCP Server 整合指南

### 10.1 什麼是 MCP？

MCP（Model Context Protocol）是由 Anthropic 於 2024 年底開源發佈的標準協定，為 AI 模型提供存取外部資料的「標準介面」。可以把 MCP 想像成 AI 的 USB-C 接口——讓 Claude 等模型直接查詢專業資料庫，而不需要手動撰寫 API 呼叫程式碼。OpenAI、Google 等也已宣佈支援 MCP。

截至 2025 年底，已有超過 15,000 個 MCP server 被部署。以下彙整本文涵蓋的各 API 對應的 MCP server，含安裝與設定方式。

### 10.2 MCP Server 總覽表

| 資料源 | MCP Server（推薦） | 語言 | 安裝方式 | API Key 需求 | GitHub |
|---|---|---|---|---|---|
| **openFDA** | OpenFDA-MCP-Server | Node.js | npm install | 建議（提升額度） | [Augmented-Nature/OpenFDA-MCP-Server](https://github.com/Augmented-Nature/OpenFDA-MCP-Server) |
| **PubMed** | pubmed-mcp-server (cyanheads) | Node.js | npx | 建議（NCBI Key） | [cyanheads/pubmed-mcp-server](https://github.com/cyanheads/pubmed-mcp-server) |
| **PubMed** | mcp-simple-pubmed | Python | pip/smithery | 建議 | [andybrandt/mcp-simple-pubmed](https://github.com/andybrandt/mcp-simple-pubmed) |
| **PubMed** | pubmed-mcp (Claude Code 最佳) | Node.js | npx (一行) | 建議 | [ncukondo/pubmed-mcp](https://github.com/ncukondo/pubmed-mcp) |
| **PubChem** | PubChem-MCP-Server (cyanheads) | Node.js | npm install | 不需要 | [cyanheads/pubchem-mcp-server](https://github.com/cyanheads/pubchem-mcp-server) |
| **PubChem** | PubChem-MCP-Server (JackKuo) | Python | pip | 不需要 | [JackKuo666/PubChem-MCP-Server](https://github.com/JackKuo666/PubChem-MCP-Server) |
| **ClinicalTrials.gov** | clinicaltrialsgov-mcp-server | Node.js | npx | 不需要 | [cyanheads/clinicaltrialsgov-mcp-server](https://github.com/cyanheads/clinicaltrialsgov-mcp-server) |
| **ClinicalTrials.gov** | ClinicalTrials-MCP-Server | Python | pip | 不需要 | [JackKuo666/ClinicalTrials-MCP-Server](https://github.com/JackKuo666/ClinicalTrials-MCP-Server) |
| **USPTO 專利** | patent_mcp_server | Python/uv | uv run | 需要（ODP Key） | [riemannzeta/patent_mcp_server](https://github.com/riemannzeta/patent_mcp_server) |
| **USPTO 商標** | trademark-mcp-server | Node.js | npx | 需要（USPTO Key） | [jordanburke/trademark-mcp-server](https://github.com/jordanburke/trademark-mcp-server) |
| **整合：BioMCP** | biomcp（PubMed + 臨床試驗 + 基因變異） | Python | uv/pip | 不需要 | [genomoncology/biomcp](https://github.com/genomoncology/biomcp) |
| **整合：Healthcare MCP** | healthcare-mcp（FDA + PubMed + 臨床試驗 + ICD-10） | Node.js | npm/dxt | 不需要 | [Cicatriiz/healthcare-mcp-public](https://github.com/Cicatriiz/healthcare-mcp-public) |
| **整合：Medical MCP** | medical-mcp（FDA + WHO + PubMed + RxNorm） | Node.js | npm -g | 不需要 | [JamesANZ/medical-mcp](https://github.com/JamesANZ/medical-mcp) |

### 10.3 Claude Code 安裝方式

#### 方法一：`claude mcp add`（推薦，最簡單）

```bash
# PubMed MCP（一行搞定）
claude mcp add pubmed-mcp --scope project -- \
  npx -y @ncukondo/pubmed-mcp --email your@email.com

# 加上 NCBI API Key（提升速率限制）
claude mcp add pubmed-mcp --scope project -- \
  npx -y @ncukondo/pubmed-mcp --email your@email.com --api-key YOUR_NCBI_KEY

# 啟用快取功能
claude mcp add pubmed-mcp --scope project -- \
  npx -y @ncukondo/pubmed-mcp --email your@email.com --cache-dir ./pubmed-cache
```

#### 方法二：編輯 `.mcp.json`（專案級）或 `~/.claude.json`（全域）

在專案根目錄建立 `.mcp.json`：

```json
{
  "mcpServers": {
    "pubmed": {
      "command": "npx",
      "args": ["-y", "@ncukondo/pubmed-mcp"],
      "env": {
        "PUBMED_EMAIL": "your@email.com",
        "PUBMED_API_KEY": "your-ncbi-api-key"
      }
    },
    "biomcp": {
      "command": "uv",
      "args": ["run", "--with", "biomcp-python", "biomcp", "run"]
    },
    "patents": {
      "command": "uv",
      "args": ["--directory", "/path/to/patent_mcp_server", "run", "patent-mcp-server"],
      "env": {
        "USPTO_API_KEY": "your-uspto-api-key"
      }
    }
  }
}
```

### 10.4 Claude Desktop 設定方式

編輯設定檔：
- **macOS**：`~/Library/Application Support/Claude/claude_desktop_config.json`
- **Windows**：`%APPDATA%\Claude\claude_desktop_config.json`

```json
{
  "mcpServers": {
    "pubmed": {
      "command": "npx",
      "args": ["-y", "@cyanheads/pubmed-mcp-server"],
      "env": {
        "NCBI_API_KEY": "your_ncbi_api_key_here"
      }
    },
    "pubchem": {
      "command": "npx",
      "args": ["-y", "@cyanheads/pubchem-mcp-server"]
    },
    "clinicaltrials": {
      "command": "npx",
      "args": ["-y", "@cyanheads/clinicaltrialsgov-mcp-server"]
    },
    "openfda": {
      "command": "node",
      "args": ["/path/to/OpenFDA-MCP-Server/build/index.js"],
      "env": {
        "FDA_API_KEY": "your_fda_api_key_here"
      }
    },
    "patents": {
      "command": "uv",
      "args": ["--directory", "/path/to/patent_mcp_server", "run", "patent-mcp-server"],
      "env": {
        "USPTO_API_KEY": "your_uspto_api_key"
      }
    },
    "trademark": {
      "command": "npx",
      "args": ["trademark-mcp-server"],
      "env": {
        "USPTO_API_KEY": "your_uspto_api_key_here"
      }
    },
    "biomcp": {
      "command": "uv",
      "args": ["run", "--with", "biomcp-python", "biomcp", "run"]
    }
  }
}
```

> ⚠️ 修改設定後須重啟 Claude Desktop 才會生效。

### 10.5 各 MCP Server 詳細說明

#### 10.5.1 openFDA MCP Server

**提供的工具（Tools）**：查詢藥物不良事件、藥品標示、召回、核准紀錄、NDC 目錄、醫材 510(k)、醫材分類、醫材不良事件等。

```bash
# 安裝
git clone https://github.com/Augmented-Nature/OpenFDA-MCP-Server.git
cd OpenFDA-MCP-Server
npm install && npm run build

# 設定 API Key（可選，提升額度）
export FDA_API_KEY="your_api_key_here"
```

使用情境：「查詢 metformin 的不良事件前 10 名反應」→ Claude 自動呼叫 openFDA 工具。

#### 10.5.2 PubMed MCP Server

**推薦方案比較**：

| Server | 特色 | 工具數 | 最適場景 |
|---|---|---|---|
| `@ncukondo/pubmed-mcp` | 一行安裝、支援快取、Claude Code 友好 | 2（search + fetch） | Claude Code 日常使用 |
| `@cyanheads/pubmed-mcp-server` | 功能最完整、支援圖表生成 | 5+ | 深度文獻研究 |
| `mcp-simple-pubmed` | Smithery 一鍵安裝、PICO 搜尋提示 | 3+ | 系統性回顧 |
| `Augmented-Nature/PubMed-MCP-Server` | 16 個工具、MeSH、作者、期刊搜尋 | 16 | 全方位文獻分析 |

**cyanheads 版提供的工具**：
- `search_pubmed` — 以關鍵字搜尋 PubMed 文獻
- `get_article_details` — 取得文章詳細資訊（摘要、作者、DOI）
- `find_related_articles` — 尋找相關文獻、引用、參考文獻
- `create_research_plan` — 建立結構化研究計畫
- `generate_pubmed_chart` — 產生統計圖表

#### 10.5.3 PubChem MCP Server

**cyanheads 版提供的工具**：
- `search_compounds` — 以名稱 / SMILES / InChIKey 搜尋化合物
- `get_compound_properties` — 取得理化性質（分子量、XLogP、TPSA 等）
- `get_compound_names` — 取得同義名、IUPAC 名稱
- `similarity_search` — 結構相似性搜尋
- `substructure_search` — 子結構搜尋
- `formula_search` — 分子式搜尋
- `get_substance_record` — 取得物質紀錄
- `get_assay_record` — 取得生物活性測定紀錄
- `get_compound_xrefs` — 跨資料庫參照（PubMed、專利等）
- `get_compound_image` — 取得 2D 結構圖

```bash
# 安裝
git clone https://github.com/cyanheads/pubchem-mcp-server.git
cd pubchem-mcp-server && npm install && npm run build
```

#### 10.5.4 ClinicalTrials.gov MCP Server

**cyanheads 版提供的工具**：
- `search_studies` — 搜尋臨床試驗（支援地理篩選）
- `get_studies` — 以 NCT ID 取得試驗詳情
- `analyze_studies` — 統計分析（最多 5,000 項試驗）
- `compare_studies` — 2–5 項試驗並排比較
- `match_patient` — 根據患者特徵匹配符合條件的試驗

```json
// Claude Desktop 設定
{
  "mcpServers": {
    "clinicaltrials": {
      "command": "npx",
      "args": ["-y", "@cyanheads/clinicaltrialsgov-mcp-server"]
    }
  }
}
```

#### 10.5.5 USPTO 專利 MCP Server

**提供 51 個工具，涵蓋 6 大資料源**：

| 資料源 | 工具功能 | 需要 API Key |
|---|---|---|
| **PPUBS**（公開搜尋） | 全文搜尋已授權專利 / 公開申請案 | 不需要 |
| **PatentsView** | 專利統計分析、發明人、受讓人查詢 | 不需要 |
| **ODP** | 申請案書目、審查歷程、文件下載 | 需要 |
| **PTAB** | IPR/PGR/申訴案裁決 | 需要 |
| **Office Actions** | 審查意見通知書 | 需要 |
| **Litigation** | 專利訴訟資料 | 需要 |

```bash
# 安裝
git clone https://github.com/riemannzeta/patent_mcp_server.git
cd patent_mcp_server

# 安裝 uv（如尚未安裝）
curl -LsSf https://astral.sh/uv/install.sh | sh

# 設定 .env
echo "USPTO_API_KEY=your_key_here" > .env

# Claude Code 整合
claude mcp add patents --scope project -- \
  uv --directory /path/to/patent_mcp_server run patent-mcp-server
```

#### 10.5.6 USPTO 商標 MCP Server

**提供的工具**：
- `trademark_search_by_serial` — 以序號搜尋商標
- `trademark_search_by_registration` — 以註冊號搜尋
- `trademark_status` — 取得商標狀態（HTML 格式）

```bash
# 直接使用 npx
export USPTO_API_KEY=your_api_key_here
npx trademark-mcp-server

# 或 Docker
docker run -d -p 3000:3000 -e USPTO_API_KEY=your_key \
  ghcr.io/jordanburke/trademark-mcp-server:latest
```

#### 10.5.7 BioMCP（整合型 — ★ 強烈推薦）

BioMCP 是一站式生醫 MCP 工具包，整合 PubMed/PubTator3、bioRxiv/medRxiv、ClinicalTrials.gov、NCI 臨床試驗、MyVariant.info、cBioPortal、OncoKB 等多個資料源，提供 **21 個工具**。

```bash
# 安裝（一行）
pip install biomcp-python
# 或
uv pip install biomcp-python

# 命令列直接使用
biomcp article search --gene BRAF --disease Melanoma
biomcp trial search --condition "Lung Cancer" --phase PHASE3
biomcp variant search --gene TP53 --significance pathogenic

# Claude Code 整合
claude mcp add biomcp --scope project -- \
  uv run --with biomcp-python biomcp run

# Claude Desktop 設定
# 加入以下至 claude_desktop_config.json：
```

```json
{
  "mcpServers": {
    "biomcp": {
      "command": "uv",
      "args": ["run", "--with", "biomcp-python", "biomcp", "run"]
    }
  }
}
```

**BioMCP 主要工具**：
| 工具 | 功能 |
|---|---|
| `article_searcher` | 搜尋 PubMed/PubTator3 及預印本 |
| `article_getter` | 取得文章詳細資訊（支援 PMID 及 DOI） |
| `trial_searcher` | 搜尋 ClinicalTrials.gov 或 NCI |
| `trial_getter` | 取得試驗完整詳情 |
| `trial_protocol_getter` | 取得試驗方案資訊 |
| `trial_outcomes_getter` | 取得試驗結果指標 |
| `trial_locations_getter` | 取得試驗執行地點 |
| `variant_searcher` | 搜尋基因變異 |
| `variant_getter` | 取得變異詳細注釋（ClinVar、COSMIC、dbSNP 等） |
| `think` | 結構化推理工具 |

#### 10.5.8 Healthcare MCP（整合型）

涵蓋 FDA、PubMed、ClinicalTrials.gov、ICD-10、medRxiv、NCBI Bookshelf，提供一鍵安裝的 DXT 擴充檔。

```bash
# 安裝
npm install -g healthcare-mcp
# 或下載 healthcare-mcp.dxt 一鍵安裝至 Claude Desktop
```

### 10.6 MCP 使用情境範例

#### 情境一：藥物安全性調查（FDA + PubMed + ClinicalTrials.gov）

> 用戶提問：「pembrolizumab 有哪些嚴重不良反應？最近有沒有相關研究？正在進行哪些新試驗？」

Claude 透過 MCP 自動執行：
1. 呼叫 **openFDA MCP** → 查詢藥物不良事件統計
2. 呼叫 **PubMed MCP** → 搜尋最近 2 年的安全性文獻
3. 呼叫 **ClinicalTrials MCP** → 列出正在招募的 Phase III 試驗

#### 情境二：專利佈局分析（USPTO Patent + PubChem）

> 用戶提問：「分析 CRISPR 基因編輯領域的美國專利趨勢，以及相關化合物的活性數據。」

Claude 透過 MCP 自動執行：
1. 呼叫 **Patent MCP** → PPUBS 搜尋 CRISPR 相關專利
2. 呼叫 **PatentsView** 工具 → 統計年度趨勢、主要受讓人
3. 呼叫 **PubChem MCP** → 查詢相關化合物的生物活性

#### 情境三：一站式生醫研究（BioMCP）

> 用戶提問：「BRAF V600E 突變在黑色素瘤中的治療現狀？」

Claude 透過 BioMCP 自動執行：
1. `think` → 結構化分析問題
2. `article_searcher` → 搜尋文獻
3. `trial_searcher` → 搜尋相關臨床試驗
4. `variant_getter` → 取得 BRAF V600E 變異注釋

### 10.7 MCP vs 直接 API 比較

| 面向 | 直接呼叫 API | 使用 MCP Server |
|---|---|---|
| **適合場景** | 批次處理、自動化腳本、大量下載 | AI 對話中即時查詢、互動式研究 |
| **設定複雜度** | 需自行處理認證、解析、錯誤處理 | 一次設定，Claude 自動使用 |
| **靈活度** | 最高（完整控制所有參數） | 受限於 MCP Server 暴露的工具 |
| **速率控制** | 需自行管理 | MCP Server 內建速率限制 |
| **Claude Code 整合** | 透過 `requests` / `curl` 手動呼叫 | 原生整合，自然語言即可觸發 |
| **建議** | 寫程式/pipeline 時使用 | 互動式研究/探索時使用 |

---

## 附錄：快速參考卡

### 「我想查⋯」→ 該用哪個 API / MCP？

| 需求 | 建議 API | 對應 MCP Server |
|---|---|---|
| 藥物核准資訊 | openFDA `/drug/drugsfda.json` | OpenFDA MCP / Healthcare MCP |
| 藥物副作用/不良事件 | openFDA `/drug/event.json` | OpenFDA MCP |
| 藥品標示/仿單 | openFDA `/drug/label.json` | OpenFDA MCP |
| 醫材 510(k) / PMA 核准 | openFDA `/device/510k.json` 或 `/device/pma.json` | OpenFDA MCP |
| 醫材不良事件 | openFDA `/device/event.json` | OpenFDA MCP |
| 食品/藥品/醫材召回 | openFDA `/{category}/enforcement.json` | OpenFDA MCP |
| 化合物結構/性質 | PubChem PUG-REST | PubChem MCP (cyanheads) |
| 化合物生物活性 | PubChem PUG-REST + `/assaysummary` | PubChem MCP (cyanheads) |
| 醫學文獻搜尋 | NCBI E-utilities (PubMed) | PubMed MCP / BioMCP |
| 全文文獻取得 | NCBI E-utilities (PMC) | PubMed MCP (Augmented-Nature) |
| 基因/蛋白質序列 | NCBI E-utilities (`db=gene`, `db=protein`) | （暫無專用 MCP） |
| 基因變異注釋 | MyVariant.info | BioMCP |
| 臨床試驗搜尋 | ClinicalTrials.gov v2 | ClinicalTrials MCP / BioMCP |
| NIH 研究經費/補助 | NIH RePORTER | （暫無專用 MCP） |
| 美國專利搜尋/分析 | USPTO PatentsView | Patent MCP Server |
| 專利申請審查進度 | USPTO ODP | Patent MCP Server |
| 專利訴訟/PTAB | USPTO PTAB + Litigation | Patent MCP Server |
| 商標狀態查詢 | USPTO TSDR | Trademark MCP Server |
| 專利大量資料 | USPTO Bulk Data / PatentsView Download | （建議直接下載） |
| **一站式生醫查詢** | 組合多個 API | **BioMCP（★ 推薦）** |
| **一站式醫療資訊** | 組合多個 API | **Healthcare MCP** |
