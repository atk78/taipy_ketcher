# Taipy Ketcher

## 開発環境

- Python: 3.12.11

ライブラリ

- flask: 3.1.1
- rdkit: 2025.3.3
- taipy: 4.1.0

本番は`gunicorn`などを用いる

## 前準備

[ketcher-standalone](https://github.com/epam/ketcher/tree/master)を取得。`standalone`ディレクトリ名を`ketcher`に変更し`src/ketcher`に設置。また、`src/ketcher/index.html`に以下のスクリプトを追記した。

```html
    <script>
    window.addEventListener('DOMContentLoaded', function() {
        function waitForKetcher() {
            if (window.ketcher && window.ketcher.getSmiles) {
                // 分子編集イベントにフック
                window.ketcher.editor?.subscribe('change', function() {
                    window.ketcher.getSmiles().then(function(smiles) {
                        // 空分子（smiles=""）は送信しない
                        if (smiles && smiles.length > 0) {
                            // 環境に応じて"http://localhost:8000/api/smiles"を変更
                            fetch("http://localhost:8000/api/smiles", {
                                method: "POST",
                                headers: { "Content-Type": "application/json" },
                                body: JSON.stringify({ smiles: smiles })
                            });
                        }
                    });
                });
            } else {
                setTimeout(waitForKetcher, 500);
            }
        }
        waitForKetcher();
    });
    </script>
```

また、Pythonスクリプト中で

```python
@app.route("/api/smiles", methods=["POST"])
def receive_smiles():
    data = request.json
    try:
        smiles = Chem.MolToSmiles(
            Chem.MolFromSmiles(data.get("smiles", ""))
            , canonical=True, kekuleSmiles=False
        )
        Data.set_smiles(smiles)
    except Exception as e:
        Data.set_smiles("")
    gui.broadcast_change("smiles", smiles)
    return jsonify({"status": "ok"})
```

によりKetcherからデータを取得する。
