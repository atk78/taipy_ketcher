from taipy.gui import Gui
from flask import Flask, request, jsonify
from flask_cors import CORS
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
import threading

# スレッドセーフな状態管理
class Data:
    _lock = threading.Lock()
    _smiles = None

    @classmethod
    def set_smiles(cls, smiles):
        with cls._lock:
            cls._smiles = smiles

    @classmethod
    def get_smiles(cls):
        with cls._lock:
            return cls._smiles

# *** Flask APIサーバー定義 ***
app = Flask(__name__)
# KetcherからのCORSリクエストを許可
CORS(app)

# POSTリクエストでSMILESを受け取るエンドポイント
@app.route("/api/smiles", methods=["POST"])
def receive_smiles():
    data = request.json
    # SMILESを受け取ってDataクラスに保存
    try:
        smiles = Chem.MolToSmiles(
            Chem.MolFromSmiles(data.get("smiles", ""))
            , canonical=True, kekuleSmiles=False
        )
        Data.set_smiles(smiles)
    except Exception as e:
        Data.set_smiles("")
    # GUIに変更を通知
    gui.broadcast_change("smiles", smiles)
    # ステータスをJSON形式で返す
    return jsonify({"status": "ok"})

# Flaskサーバーを起動する関数
def run_flask():
    app.run(port=8000)

# *** Taipy GUI側 ***
smiles = None

def on_init(state):
    Data.set_smiles("")
    state.smiles = Data.get_smiles()

def refresh_smiles(state):
    Data.set_smiles(state.smiles)

def calc_property(state):
    # ここでSMILESからプロパティを計算する処理を実装
    # 今回はダミーとして、SMILESの長さを返す
    mol = Chem.MolFromSmiles(state.smiles)
    if mol is not None:
        molwt = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)  # ダミー分子のLogP
        print(f"SMILES: {state.smiles}, 分子量: {molwt}, LogP: {logp}")

# *** ページ定義 ***

page = """
# Taipy Ketcher

<iframe src="/src/ketcher/index.html" width="900" height="700"></iframe>

**SMILES:** <|{smiles}|>

<|ボタン|button|on_action=calc_property|label=「計算」|>
"""


if __name__ == "__main__":
    # Flaskサーバーを別スレッドで立ち上げる
    flask_thread = threading.Thread(target=run_flask, daemon=True)
    flask_thread.start()

    gui = Gui(page=page)
    gui.run(
        title="Taipy Ketcher",
        port=8080,
        run_browser=False,  # ブラウザは自動で開かないようにする
        on_init=on_init,
        static_folder="src/ketcher",  # Ketcherの静的ファイル
        static_url_path="/src/ketcher"  # URLパスを指定
    )
