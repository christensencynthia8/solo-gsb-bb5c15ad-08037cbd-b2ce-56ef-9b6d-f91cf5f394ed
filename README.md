# DNA 限制性酶切模拟 API

基于 **Python + FastAPI + SQLite** 的限制性内切酶消化模拟服务。提交线性或环形 DNA
序列（仅 IUPAC 碱基）以及多种酶的识别模式与正/反链切口偏移，服务在**两个方向**
查找识别位点（允许重叠、含简并碱基、可跨越环形序列首尾），对同一物理切口去重后
按拓扑结构生成片段，返回坐标、长度、两端的平端/黏端类型及突出端序列。

## 快速启动

一条指令完成构建、依赖安装与数据库初始化：

```bash
docker compose up --build
```

依赖在镜像构建时安装（`pip install -r requirements.txt`），SQLite 数据库在容器
启动时自动建库建表，无需任何手工步骤。启动后访问 **http://localhost:8000/docs**
即可看到包含全部接口的交互式文档。

其它常用 Docker 操作：

```bash
# 后台运行
docker compose up --build -d

# 查看日志 / 停止
docker compose logs -f
docker compose down

# 仅构建镜像
docker build -t digest-api .
docker run --rm -p 8000:8000 -v "$PWD/data:/app/data" digest-api

# 在容器内运行测试套件
docker build -t digest-api .
docker run --rm digest-api sh -c "pip install --no-cache-dir pytest httpx && pytest"
```

## 配置

| 配置项 | 默认值 | 说明 |
| --- | --- | --- |
| 环境变量 `DATABASE_PATH` | `/app/data/digest.db`（容器内） | SQLite 数据库文件路径，目录不存在时自动创建 |
| 端口映射 | `8000:8000`（compose 中 `ports`） | 宿主机端口：容器端口 |
| 数据卷 | `digest-data` → `/app/data` | 持久化数据库文件 |

不使用 Docker 时本地运行：

```bash
pip install -r requirements.txt
uvicorn app.main:app --host 0.0.0.0 --port 8000
```

## API 地址

- 服务地址：**http://localhost:8000**
- 交互式文档（Swagger UI）：**http://localhost:8000/docs**
- 健康检查：`GET /health`

### 接口一览

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| `POST` | `/digests` | 提交序列与酶列表，执行酶切、保存并返回完整结果（201） |
| `GET` | `/digests` | 列出全部提交的摘要（按 id 升序） |
| `GET` | `/digests/{id}` | 按 id 查询保存的完整酶切结果（不存在返回 404） |
| `POST` | `/digests/{id}/variant-analysis` | 对已保存记录施加一组单碱基替换并重跑酶切，返回切口与片段长度差异（不落库，记录不存在返回 404） |
| `POST` | `/digests/{id}/enzyme-panel-plan` | 枚举候选酶组合，寻找达到目标片段数的最少酶集合（只读不落库，记录不存在返回 404，无解返回 422） |
| `POST` | `/digests/{id}/partial-digest-distribution` | 为每个已保存切口指定独立切割概率，枚举切/不切组合并按片段长度多重集合并概率分布（只读不落库，记录不存在返回 404，坐标不匹配或切口超过 12 个返回 422） |
| `GET` | `/health` | 健康检查 |

### 提交示例

```bash
curl -X POST http://localhost:8000/digests \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "AAAAGAATTCTTTT",
    "topology": "linear",
    "enzymes": [
      {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5},
      {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}
    ]
  }'
```

- `sequence`：仅允许 IUPAC 碱基 `ACGTRYSWKMBDHVN`（大小写均可，存储为大写）。
- `topology`：`linear`（线性）或 `circular`（环形）。
- `enzymes`：至少一种酶。
  - `pattern`：识别模式，同样仅允许 IUPAC 碱基；
  - `cut_top` / `cut_bottom`：正/反链切口偏移，以识别位点 5' 端为基准，
    合法范围为 **`0 ~ len(pattern)`**。

响应（节选）：

```json
{
  "id": 1,
  "topology": "linear",
  "sequence_length": 14,
  "cuts": [
    {
      "top_cut": 5,
      "bottom_cut": 9,
      "end_type": "5prime",
      "overhang": "AATT",
      "enzymes": ["EcoRI"],
      "sites": [{"enzyme": "EcoRI", "start": 4, "strand": "+", "matched": "GAATTC"}]
    }
  ],
  "fragments": [
    {"index": 0, "start": 0, "end": 5, "length": 5, "sequence": "AAAAG",
     "wraps_origin": false,
     "left_end": {"type": "none", "overhang": ""},
     "right_end": {"type": "5prime", "overhang": "AATT"}},
    {"index": 1, "start": 5, "end": 14, "length": 9, "sequence": "AATTCTTTT",
     "wraps_origin": false,
     "left_end": {"type": "5prime", "overhang": "AATT"},
     "right_end": {"type": "none", "overhang": ""}}
  ]
}
```

## 变异影响分析

对已保存的酶切记录施加一组 **0-based 单碱基替换**，在**内存中**用变异后序列
重跑酶切（拓扑与酶列表沿用原记录），并与原结果逐项对比。**分析过程不写库**：
原记录保持不变，也不会产生新的提交。

### 请求

```bash
curl -X POST http://localhost:8000/digests/1/variant-analysis \
  -H 'Content-Type: application/json' \
  -d '{"variants": [{"position": 4, "alt": "C"}]}'
```

- `variants`：至少一个替换；每项 `position` 为 0-based 位置，`alt` 为单个
  IUPAC 碱基（大小写均可，存储为大写）。
- 允许同一碱基的"同义替换"（`alt` 与原碱基相同），此时差异全部为空。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| `position` 越界（`>= 序列长度`，含上界） | **422** |
| `position` 为负数 / `alt` 为非 IUPAC 碱基 | **422** |
| 同一 `position` 重复出现 | **422** |
| `variants` 为空列表 | **422** |

### 响应与差异口径

```json
{
  "id": 1,
  "sequence_length": 14,
  "variant_sequence": "AAAACAATTCTTTT",
  "variants": [{"position": 4, "ref": "G", "alt": "C"}],
  "cuts": {
    "added": [],
    "removed": [ {"top_cut": 5, "bottom_cut": 9, "...": "（完整切口对象）"} ],
    "changed": []
  },
  "fragments": {
    "before_lengths": [5, 9],
    "after_lengths": [14],
    "added": [],
    "removed": [ {"start": 5, "...": "（完整片段对象）"} ],
    "changed": [
      {"start": 0, "before_length": 5, "after_length": 14, "length_delta": 9}
    ]
  }
}
```

- **切口身份 = 物理坐标 `(top_cut, bottom_cut)`**：
  - `added`：仅变异后存在的切口（**新增**）；
  - `removed`：仅原结果存在的切口（**消失**）；
  - `changed`：坐标相同、但末端或归属发生变化，`before` / `after` 为两个完整
    切口对象，`changes` 列出变化字段（`end_type`、`overhang`、`enzymes`、
    `sites` 的子集）——例如简并位点存活但突出端序列被替换改变；
  - 切口"移动位置"（坐标/末端位置改变）表现为一条 `removed` 加一条 `added`。
- **片段差异**：替换不改变序列长度，片段按骨架起点 `start` 配对。`changed`
  给出起点相同片段的长度差异（`length_delta = after_length - before_length`）；
  `added` / `removed` 为仅一侧存在的完整片段；`before_lengths` /
  `after_lengths` 为物理顺序的片段长度清单，便于快速对比。
- **稳定排序**：回显变异按 `position` 升序（与传入顺序无关）；切口差异按
  `(top_cut, bottom_cut)` 排序；片段差异按 `start` 排序，相同请求结果必然一致。

## 酶组合面板规划

基于已保存记录的**序列与拓扑**，枚举请求中给出的候选酶（**至多 8 种**）的所有
组合并复用酶切计算，寻找**恰好达到目标片段数的最少酶集合**。**规划过程不写库**：
原记录保持不变，也不会产生新的提交。

### 请求

以下示例假设 `id=1` 的记录序列为 `AAAAGAATTCTTTTGGATCCGGGG`（线性，含 EcoRI、
BamHI 各一个位点；可先通过 `POST /digests` 创建）：

```bash
curl -X POST http://localhost:8000/digests/1/enzyme-panel-plan \
  -H 'Content-Type: application/json' \
  -d '{
    "target_fragment_count": 3,
    "candidates": [
      {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5},
      {"name": "BamHI", "pattern": "GGATCC", "cut_top": 1, "cut_bottom": 5},
      {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}
    ]
  }'
```

- `target_fragment_count`：目标片段数，整数且 `>= 1`。
- `candidates`：1~8 种候选酶，字段与校验规则和提交酶切时的 `enzymes` 完全一致；
  完全相同（名称/模式/偏移均相同）的候选会自动去重。

### 响应

```json
{
  "id": 1,
  "target_fragment_count": 3,
  "selected_enzymes": [
    {"name": "BamHI", "pattern": "GGATCC", "cut_top": 1, "cut_bottom": 5},
    {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
  ],
  "fragment_count": 3,
  "cuts": [
    {"top_cut": 5, "bottom_cut": 9, "end_type": "5prime", "overhang": "AATT",
     "enzymes": ["EcoRI"], "sites": ["..."]},
    {"top_cut": 15, "bottom_cut": 19, "end_type": "5prime", "overhang": "GATC",
     "enzymes": ["BamHI"], "sites": ["..."]}
  ]
}
```

- `selected_enzymes`：选中的最少酶集合，按酶名排序返回；`fragment_count` 为其实际
  达到的片段数（恒等于目标值）；`cuts` 为该组合切口的完整对象（含坐标、末端类型、
  归属酶与位点），与 `POST /digests` 结果中的切口结构一致。

### 择优规则与状态码

- **最少酶优先**：按组合大小从小到大枚举，第一个达到目标片段数的规模即为答案；
  空集也参与枚举——不切割时线性/环形分子均为 1 条片段，因此目标为 1 时可能返回
  空的 `selected_enzymes`。
- **酶数相同时按酶名序列择优**：比较排序后的酶名元组，字典序最小者胜出；候选在
  请求中的传入顺序不影响结果，相同请求必然得到相同方案。
- **无解返回 422**：所有候选组合（含空集）都无法恰好达到目标片段数时，返回
  `422` 并说明原因。

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| 所有候选组合均无法达到目标片段数 | **422** |
| `target_fragment_count` < 1 | **422** |
| `candidates` 为空或超过 8 种 | **422** |
| 候选酶非法（非 IUPAC 模式 / 切口偏移越界） | **422** |

## 部分酶切概率分布

部分酶切（partial digest）模拟：为已保存记录中的**每个物理切口**指定
`0 ~ 1` 的独立切割概率，枚举全部"切开/未切开"组合（**至多 12 个切口**，
即 4096 种），复用线性/环形片段重建逻辑计算每种组合的片段长度多重集；
多重集相同的结果**合并并累加概率**，按概率降序返回完整分布。
**分析过程不写库**：原记录保持不变，也不会产生新的提交。

服务启动后访问地址为 **http://localhost:8000**，交互式文档见
**http://localhost:8000/docs**（其中可在线试调本接口）。

### 请求

以下示例假设 `id=1` 的记录为线性序列 `AAAAGAATTCTTTT`、含一个 EcoRI
切口 `(top_cut=5, bottom_cut=9)`（可先通过 `POST /digests` 创建）：

```bash
curl -X POST http://localhost:8000/digests/1/partial-digest-distribution \
  -H 'Content-Type: application/json' \
  -d '{"cuts": [{"top_cut": 5, "bottom_cut": 9, "probability": 0.75}]}'
```

- `cuts`：为**每个**已保存物理切口指定概率，切口以物理坐标
  `(top_cut, bottom_cut)` 标识（与 `GET /digests/{id}` 返回的切口一一对应，
  传入顺序不限），`probability` 取值 `0 ~ 1`。
- 坐标集合必须与已保存切口**完全一致**：缺失、多余或重复均返回 422。

### 响应

```json
{
  "id": 1,
  "topology": "linear",
  "cut_count": 1,
  "outcome_count": 2,
  "total_probability": 1.0,
  "outcomes": [
    {"fragment_lengths": [5, 9], "fragment_count": 2, "probability": 0.75},
    {"fragment_lengths": [14], "fragment_count": 1, "probability": 0.25}
  ]
}
```

- `outcomes`：合并后的分布，按 `probability` **降序**排列（概率相同按片段
  长度多重集升序，输出完全确定）；`fragment_lengths` 为升序的片段长度多重集。
- **合并口径**：只要片段长度多重集相同即合并——例如线性分子上只切 5 位得
  `[5, 9]` 与只切 9 位得 `[9, 5]`；环形分子 0 个与 1 个切口均为全长片段。
- `total_probability`：全部结果概率之和，与 1 的误差不超过 `1e-9`；
  概率为 0 的组合（某切口概率取 0 或 1 时）不列入分布。

### 校验与状态码

| 情形 | 状态码 |
| --- | --- |
| 记录不存在 | **404** |
| 已保存切口超过 12 个（枚举上限） | **422** |
| 坐标缺失 / 多余（与已保存切口不一致） | **422** |
| 坐标重复 | **422** |
| `probability` 超出 `[0, 1]` | **422** |

## 结果字段约定

- **坐标**：均为 0-based，作用于 top 链（提交序列方向）；切口位置 `c` 表示断裂
  位于第 `c-1` 与第 `c` 个碱基之间。片段区间为半开区间 `[start, end)`。
- **末端类型 `end_type` / `left_end.type` / `right_end.type`**：
  - `blunt`：平端（正/反链切口位置相同）；
  - `5prime`：5' 黏端；`3prime`：3' 黏端；
  - `none`：线性分子的天然末端，或未切割的完整分子。
- **`overhang`**：突出端序列，按 top 链 5'→3' 方向取自两切口之间的区段；
  平端与 `none` 时为空字符串。
- **环形序列**：切口位置对序列长度取模；`wraps_origin: true` 表示该片段跨越
  首尾原点，其 `end` 为取模后的坐标（可能小于 `start`）。单个切口开环时返回
  一条全长片段；无切口时保留整条环状分子。
- **结果顺序**：切口按 `(top_cut, bottom_cut)` 排序，片段按起点排序，酶名与
  位点列表均排序——**相同输入的结果顺序必然一致**。

## 行为规则

- **双向搜索**：同时在 top 链搜索 `pattern` 与 `reverse_complement(pattern)`
  （即酶结合于反链的情形）；回文模式只搜索一次。
- **重叠位点**：允许识别位点相互重叠，每个命中都会计算切口。
- **简并碱基**：模式与序列中的 IUPAC 简并碱基按"可能集合有交集"匹配
  （如模式中的 `A` 可匹配序列中的 `R`）。
- **跨原点**：环形序列在扩展序列上搜索，跨越首尾的识别位点与切口均可命中。
- **物理切口去重**：不同酶/不同位点若产生相同的 `(top_cut, bottom_cut)`，
  只保留一条切口记录，酶名合并进 `enzymes` 列表。
- **非法偏移整体拒绝**：任一酶的切口偏移超出 `[0, len(pattern)]`（含负数），
  或序列/模式含非 IUPAC 字符，整个请求以 **422** 拒绝，**不会有任何数据落库**。
- **无切口**：返回一条覆盖整条分子的片段，两端类型为 `none`。

## 项目结构

```
app/
├── main.py      # FastAPI 入口与路由
├── digest.py    # 酶切核心算法（搜索、去重、片段生成）
├── variant.py   # 变异施加与切口/片段差异对比
├── panel.py     # 酶组合面板规划（枚举候选组合、最少酶择优）
├── partial.py   # 部分酶切概率分布（枚举切/不切组合、按长度多重集合并）
├── iupac.py     # IUPAC 编码表与反向互补
├── schemas.py   # 请求/响应模型与校验
└── db.py        # SQLite 持久化（启动时自动建表）
tests/
├── test_digest.py       # 酶切算法测试
├── test_api.py          # 提交/查询/校验 API 测试
├── test_variant_api.py  # 变异影响分析 API 测试
├── test_panel_api.py    # 酶组合面板规划 API 测试
└── test_partial_api.py  # 部分酶切概率分布 API 测试
Dockerfile
docker-compose.yml
requirements.txt
requirements-dev.txt
```

## 运行测试

本地运行：

```bash
pip install -r requirements-dev.txt
pytest
```

容器内运行（无需本地 Python 环境）：

```bash
docker build -t digest-api .
docker run --rm digest-api sh -c "pip install --no-cache-dir -r requirements-dev.txt && pytest"
```
