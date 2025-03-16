#!/bin/bash
# docking_automationパッケージのテストを実行するスクリプト

# スクリプトの場所を取得
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

cd "$PROJECT_ROOT"

# 実行モード
MODE=${1:-"all"}  # 引数がなければ "all" をデフォルトとする

echo "==============================================="
echo "テスト実行モード: $MODE"
echo "プロジェクトルート: $PROJECT_ROOT"
echo "==============================================="

# 環境変数
export PYTHONPATH="$PROJECT_ROOT:$PYTHONPATH"
export MOCK_DOCKING=1  # モックモードを有効化

# テストの実行
case $MODE in
  "unit")
    echo "ユニットテストを実行します..."
    pytest -xvs docking_automation/tests/unit/
    ;;
  "integration")
    echo "統合テストを実行します..."
    pytest -xvs docking_automation/tests/integration/
    ;;
  "scenario")
    echo "シナリオテストを実行します..."
    pytest -xvs docking_automation/tests/scenario/
    ;;
  "all" | *)
    echo "すべてのテストを実行します..."
    pytest -xvs docking_automation/tests/
    ;;
esac

exit $?