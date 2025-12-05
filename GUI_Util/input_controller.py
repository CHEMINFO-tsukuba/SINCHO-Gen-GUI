import streamlit as st
import subprocess
import os
import tempfile
import io
from Bio.PDB import PDBParser
from fractions import Fraction
import py3Dmol
from streamlit.components.v1 import html
import yaml
import shutil
import glob
import pandas as pd
import numpy as np



class InputController:
    def __init__(self):
        st.write("")

    
    def process(self, sub_tab):


        #"General", "Upload Complex", "Select Hit Ligand", "MD Settings", "SINCHO Settings", "ChemTS Settings", "AAScore Settings", "Summary"

        if sub_tab == "General":
            st.title("General Settings")

            st.session_state.general_settings.setdefault("dir_step", "init")

            st.session_state.general_settings["use_num_threads"] = st.number_input(
                "並列スレッド数", value=12, step=1
            )
            st.session_state.general_settings["directory"] = st.text_input(
                "出力ディレクトリ", value="./output"
            )

            wdir = os.path.join(os.getcwd(), st.session_state.general_settings["directory"])

            if st.session_state.general_settings["dir_step"] == "done":
                st.success(f"ディレクトリが設定されました: {wdir}")
                # 次の設定UIなどをここに書ける
                st.write("ここからさらに操作を続けられます")

            elif st.session_state.general_settings["dir_step"] == "init":
                if st.button("上記のディレクトリに設定する"):
                    if os.path.exists(wdir):
                        st.session_state.general_settings["dir_step"] = "confirm"
                        try:
                            st.rerun()
                        except Exception as e:
                            st.experimental_rerun()
                    else:
                        try:
                            os.makedirs(os.path.join(wdir, "99_TMP"), exist_ok=True)
                            st.session_state.general_settings["tmp_dir"] = os.path.join(wdir, "99_TMP")
                            st.session_state.general_settings["dir_step"] = "done"
                            try:
                                st.rerun()
                            except Exception as e:
                                st.experimental_rerun()
                        except Exception as e:
                            st.error(f"作成に失敗しました: {e}")

            elif st.session_state.general_settings["dir_step"] == "confirm":
                st.warning(f"ディレクトリ {wdir} は既に存在します。")
                if st.button("上書きを許可する"):
                    try:
                        os.makedirs(os.path.join(wdir, "99_TMP"), exist_ok=True)
                        st.session_state.general_settings["tmp_dir"] = os.path.join(wdir, "99_TMP")
                        st.session_state.general_settings["dir_step"] = "done"
                        try:
                            st.rerun()
                        except Exception as e:
                            st.experimental_rerun()
                    except Exception as e:
                        st.error(f"作成に失敗しました: {e}")

        if sub_tab == "Initial Upload":
            st.title("Initial Complex Structure Upload")


            uploaded_file = st.file_uploader("複合体のPDBファイルをアップロードしてください", type=["pdb"])

            if uploaded_file:
                # session_state に記録
                st.session_state.uploaded_pdb_file = uploaded_file
                # 3D可視化
                self._pdb_3dview(uploaded_file)
                
                # 一時フォルダに保存
                tmp_path = os.path.join(st.session_state.general_settings["tmp_dir"], uploaded_file.name)
                with open(tmp_path, "wb") as f:
                    f.write(uploaded_file.getvalue())
                st.session_state.uploaded_pdb_file.path = tmp_path
                
                st.success(f"ファイルがアップロードされました ({tmp_path})")
                st.write("ファイルを変更したい場合は再度アップロードしてください")
                """
                if st.button("アップロードしたファイルをキャンセルする"):
                    st.session_state.uploaded_pdb_file = None
                    uploaded_file = None
                    try:
                        st.rerun()
                    except Exception as e:
                        st.experimental_rerun()
                """
            else:
                if st.session_state.uploaded_pdb_file:
                    # session_state に残っているものを表示
                    self._pdb_3dview(st.session_state.uploaded_pdb_file)
                    st.success(f"既にファイル({st.session_state.uploaded_pdb_file.name})がアップロードされています。")

        if sub_tab == "Hit Residue Selection":
            st.title("Hit Residue Selection")

            try:
                # まずPDBがアップロードされているかを確認
                if not st.session_state.uploaded_pdb_file:
                    st.warning("PDB ファイルがアップロードされていません。Initial Uploadタブに戻ってください。")
                    st.stop()
                
                st.session_state.residues_list = self._residue_parser(st.session_state.uploaded_pdb_file)
                
                # session_state で選択値を保持
                if "hit_residue" not in st.session_state or st.session_state.hit_residue not in st.session_state.residues_list:
                    st.session_state.hit_residue = st.session_state.residues_list[0]

                
                selected_residue = st.selectbox(
                    "ヒット化合物の残基を選択してください(1残基のみ選択可能)",
                    st.session_state.residues_list,
                    index=st.session_state.residues_list.index(st.session_state.hit_residue)
                )
                
                # session_stateに更新
                st.session_state.hit_residue = selected_residue
                
                st.success(f"選択された残基: {st.session_state.hit_residue}")
                self._pdb_3dview(st.session_state.uploaded_pdb_file, zoomres=st.session_state.hit_residue)
                st.success("ヒット化合物が選択されました。次のステップに進んでください。")

            except Exception as e:
                st.error(f"ファイルの読み込み中にエラーが発生しました: {e}")

        if sub_tab == "PDB File Editor":
            st.title("PDB File Direct Editor")
            self._pdb_editable_board()
            
        if sub_tab == "MD Settings":
            st.title("Molecular Dynamics Simulation Settings")

            if "md_settings" not in st.session_state:
                st.session_state.md_settings = {}
            if "force_field" not in st.session_state.md_settings or len(st.session_state.md_settings["force_field"]) != 4:
                st.session_state.md_settings["force_field"] = ["ff14SB", "gaff2", "tip3p", "OL3"]
            if "temperature" not in st.session_state.md_settings:
                st.session_state.md_settings["temperature"] = 300
            if "additional_parameters" not in st.session_state.md_settings:
                st.session_state.md_settings["additional_parameters"] = None
            if "box_shape" not in st.session_state.md_settings:
                st.session_state.md_settings["box_shape"] = "rectangular"
            if "box_size" not in st.session_state.md_settings:
                st.session_state.md_settings["box_size"] = 75.0
            if "buffer" not in st.session_state.md_settings:
                st.session_state.md_settings["buffer"] = 15.0
            if "snapshots" not in st.session_state.md_settings:
                st.session_state.md_settings["snapshots"] = None
            if "pr_run_time" not in st.session_state.md_settings:
                st.session_state.md_settings["pr_run_time"] = 10
            if "pr_rec_interval" not in st.session_state.md_settings:
                st.session_state.md_settings["pr_rec_interval"] = 2

            with st.expander("力場パラメータ設定"):
                st.session_state.md_settings["force_field"] = [st.selectbox("タンパク質",["ff14SB", "ff99SB", "ff19SB"], index=["ff14SB", "ff99SB", "ff19SB"].index(st.session_state.md_settings["force_field"][0])),
                                                               st.selectbox("化合物",["gaff2", "gaff"], index=["gaff2", "gaff"].index(st.session_state.md_settings["force_field"][1])),
                                                               st.selectbox("水分子",["tip3p", "spc"], index=["tip3p", "spc"].index(st.session_state.md_settings["force_field"][2])),
                                                               st.selectbox("RNA (if any)",["OL3", "OL4"], index=["OL3", "OL4"].index(st.session_state.md_settings["force_field"][3]))]
            st.success(f"選択された力場: {st.session_state.md_settings['force_field']}")
            with st.expander("シミュレーション条件"): # （平衡化過程等は固定値を使用します。今後軽量版平衡化も選択可能にする予定。）
                st.session_state.md_settings["temperature"] = st.number_input("温度 (K)(動的変数未実装：現状300K固定です)", value=st.session_state.md_settings["temperature"], step=1)
                
                st.write("化合物のパラメータファイルの追加アップロード（残基ごと）➡無い場合はGasteiger chargeを使用")
                selected_residues = st.multiselect(        
                    "パラメータを設定したい残基を選んでください（複数可）",
                    st.session_state.residues_list,
                    default=[st.session_state.hit_residue])

                # 残基 → ファイルリストの辞書を初期化
                if "additional_parameters" not in st.session_state.md_settings:
                    st.session_state.md_settings["additional_parameters"] = {}
                tmp_dir = st.session_state.general_settings["tmp_dir"]

                # 各残基についてアップロード UI と保存処理
                saved_paths = []
                for resname in selected_residues:
                    st.markdown(f"### 🔹 残基 `{resname}` のパラメータファイル")
                    files = st.file_uploader(
                        f"{resname} に対応する .prep / .frcmod ファイルをアップロード",
                        type=["prep", "frcmod"],
                        accept_multiple_files=True,
                        key=f"uploader_{resname}"
                    )

                    if files:
                        for file in files:
                            tmp_path = os.path.join(tmp_dir, resname.split(" ")[0]+os.path.splitext(file.name)[1])
                            with open(tmp_path, "wb") as f:
                                f.write(file.getvalue())
                            saved_paths.append(resname.split(" ")[0]+os.path.splitext(file.name)[1])

                # 保存
                st.session_state.md_settings["additional_parameters"] = saved_paths

                st.success(f"{selected_residues} に対応するファイルを保存しました")
                st.success(glob.glob(os.path.join(tmp_dir, "*.prep")) + glob.glob(os.path.join(tmp_dir, "*.frcmod")))




                st.session_state.md_settings["box_shape"] = st.selectbox("ボックス形状", ["rectangular", "cube"], index=["rectangular", "cube"].index(st.session_state.md_settings["box_shape"]))
                if st.session_state.md_settings["box_shape"] == "cube":
                    st.session_state.md_settings["box_size"] = st.number_input("ボックスサイズ (Å)", value=st.session_state.md_settings["box_size"], step=1.0)
                else:
                    st.session_state.md_settings["buffer"] = st.number_input("バッファサイズ (Å)", value=st.session_state.md_settings["buffer"], step=0.1)    

            st.write("Production Runの詳細設定をしてください")
            st.session_state.md_settings["pr_run_time"] = st.number_input("Production Run時間 (ns)", value=st.session_state.md_settings["pr_run_time"], step=1)
            st.session_state.md_settings["pr_rec_interval"] = st.number_input("インターバル (ns)", value= st.session_state.md_settings["pr_rec_interval"], step=1)
            snaps = Fraction(st.session_state.md_settings["pr_run_time"]) / Fraction(st.session_state.md_settings["pr_rec_interval"])
            maximum_steps = st.session_state.md_settings["pr_run_time"]/(0.002/1000000)
            if snaps.denominator==1:
                if int(snaps) > maximum_steps or int(snaps)<1:
                    st.warning(f"1 <= (Production Run時間)/(インターバル) <= {int(maximum_steps)}を満たしてください")
                else:
                    st.success(f"{snaps}コ+1コ(0ns)={snaps+1}コのスナップショットが保存され、以下の処理に使用されます。\n次のステップに進んでください。")
                    st.session_state.md_settings["snapshots"] = int(snaps)
            else:
                st.warning(f"(Production Run時間)/(インターバル)を整数値にしてください")
                st.session_state.md_settings["snapshots"] = None

        if sub_tab == "SINCHO Settings":
            st.title("SINCHO Settings")
            if "p2c_sincho_settings" not in st.session_state:
                st.session_state.p2c_sincho_settings = {}
            if "distance_range" not in st.session_state.p2c_sincho_settings:
                st.session_state.p2c_sincho_settings["distance_range"] = 10.0  # Å
            if "npairs_per_snap" not in st.session_state.p2c_sincho_settings:
                st.session_state.p2c_sincho_settings["npairs_per_snap"] = 10  # ペア数
            if "for_chemts" not in st.session_state.p2c_sincho_settings:
                st.session_state.p2c_sincho_settings["for_chemts"] = 2  # ChemTSに渡す候補数
            if "r_point_atoms" not in st.session_state.p2c_sincho_settings:
                st.session_state.p2c_sincho_settings["r_point_atoms"] = None

            if st.session_state.md_settings["snapshots"]:
                st.success(f"MDシミュレーションから得られるトラジェクトリから、\n{str(st.session_state.md_settings['snapshots']+1)}個のスナップショットを保存して以降の処理を進めます。")
                st.session_state.p2c_sincho_settings["distance_range"] = st.number_input("P2Cのポケット探索範囲(化合物からX[Å]以内のポケットのみ探索)", value=st.session_state.p2c_sincho_settings["distance_range"], step=0.1)
                st.session_state.p2c_sincho_settings["npairs_per_snap"] = st.number_input("SINCHO結果の出力数(per snapshot)", value=st.session_state.p2c_sincho_settings["npairs_per_snap"], step=1)
                st.session_state.p2c_sincho_settings["for_chemts"] = st.number_input("化合物生成に使用するSINCHO結果のペア数(per snapshot)", value=st.session_state.p2c_sincho_settings["for_chemts"], step=1)
                if st.session_state.p2c_sincho_settings["npairs_per_snap"] < st.session_state.p2c_sincho_settings["for_chemts"]:
                    st.warning("SINCHOの予測ペア数 > ChemTSに渡すペア数 である必要があります。\n設定を見直してください")
                else:
                    st.success("SINCHOの設定が完了しました。上記の設定で問題無ければ、次のステップに進んでください。")
            else:
                st.warning("スナップショット数の設定が適切ではありません。Production Run設定に戻ってください")

            st.write("R-pointsの制限（optional）")
            name_list = self._pdb_3dview_res(st.session_state.hit_residue)

            # ① 初期化：p2c_sincho_settings 側
            if st.session_state.p2c_sincho_settings["r_point_atoms"] is None:
                # 初回は「全部選択」状態にして保存
                st.session_state.p2c_sincho_settings["r_point_atoms"] = name_list

            # ② 初期化：multiselect 用の session_state
            if "selected_r_points" not in st.session_state:
                # r_point_atoms 文字列 → リスト
                r_points_str = st.session_state.p2c_sincho_settings["r_point_atoms"] or ""
                if r_points_str:
                    default_list = [x for x in r_points_str if x in name_list]
                else:
                    default_list = name_list  # 何もなければ全選択とかお好みで
                st.session_state.selected_r_points = default_list

            # ③ multiselect（ここでは default は渡さないのが重要）
            st.session_state.p2c_sincho_settings["r_point_atoms"] = st.multiselect(
                "R-pointとして使用するatomname「のみ」残してください",
                name_list,
                key="selected_r_points",
            )



        if sub_tab == "ChemTS Settings":
            st.title("ChemTSv2 Settings")
            if "chemts_settings" not in st.session_state:
                st.session_state.chemts_settings = {}
            if "num_chemts_loops" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["num_chemts_loops"] = 4  # 生成の反復回数
            if "c_val" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["c_val"] = 1.0           # C値
            if "threshold_type" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["threshold_type"] = "time"  # 終了条件のタイプ
            if "threshold" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["threshold"] = 0.05      # 終了条件の値 (時間 or 生成数)
            if "function_format" not in st.session_state.chemts_settings:
                st.session_state.chemts_settings["function_format"] = "only_sincho"  # 報酬の形式

            st.success(f"スナップショット:{st.session_state.md_settings['snapshots']+1}コ × リード展開方針:{st.session_state.p2c_sincho_settings['for_chemts']}コを用いてリード分子生成を行います。")
            with st.expander("Basic setting"):
                st.session_state.chemts_settings["num_chemts_loops"] = st.number_input("生成の反復回数", value=st.session_state.chemts_settings["num_chemts_loops"], step=1)
                st.session_state.chemts_settings["c_val"] = st.number_input("C値", value=st.session_state.chemts_settings["c_val"], step=0.1)
                st.session_state.chemts_settings["threshold_type"] = st.selectbox("1回の生成の終了条件", options=["time","generation_num"], index=["time","generation_num"].index(st.session_state.chemts_settings["threshold_type"]))
                if st.session_state.chemts_settings["threshold_type"] == "generation_num":
                    st.session_state.chemts_settings["threshold"] = st.number_input("生成数", value=st.session_state.chemts_settings["threshold"], step=1)
                elif st.session_state.chemts_settings["threshold_type"] == "time":
                    st.session_state.chemts_settings["threshold"] = st.number_input("時間 (hour)", value=st.session_state.chemts_settings["threshold"], step=0.01)
                    scaler = st.session_state.chemts_settings["threshold"]*st.session_state.chemts_settings["num_chemts_loops"]*st.session_state.p2c_sincho_settings["npairs_per_snap"]*st.session_state.md_settings["snapshots"]
                    st.write(f"ChemTSのおおよその実行時間: {round(scaler,3)}時間")
                
            with st.expander("Advanced setting"):
                st.write("後々実装予定。それまでは最後の直編集パネルでの設定をお願いします。")
            
            with st.expander("Filter setting"):
                st.write("後々実装予定。それまでは最後の直編集パネルでの設定をお願いします。")
                """
                filter_dict = {"lipinski":[["rule_of_5","rule_of_3"],0],"radical":None,"pubchem":None,"sascore":[3.5],"ring_size":[6],"pains":[["[pains_a]"],0],"donor_acceptor":None}
                for fil,op in filter_dict.items():
                    st.session_state.chemts_settings[f"use_{fil}_filter"] = st.selectbox(f"use_{fil}_filter?", options=["True","False"], index=["True","False"].index(st.session_state.chemts_settings.get(f"use_{fil}_filter")))
                    if st.session_state.chemts_settings[f"use_{fil}_filter"] == "True":
                        if op==None:
                            pass
                        elif len(op)==1:
                            st.session_state.chemts_settings[f"{fil}_threshold"] = st.number_input(f"____{fil}_threshold",value=st.session_state.chemts_settings[f"{fil}_threshold"] ,step=0.1)
                        elif len(op)==2:
                            st.session_state.chemts_settings[f"{fil}_type"] = st.selectbox(f"____{fil}_type", options=op[0], index=op[1])
                """

            with st.expander("Reward setting"):
                st.session_state.chemts_settings["function_format"] = st.selectbox("報酬の形式を選択してください", options=["only_sincho", "cns"], index=["only_sincho", "cns"].index(st.session_state.chemts_settings["function_format"]))
                st.write("各項の形式は固定値を使用します。修正は最後の直編集パネルで行ってください")
                st.write("現状、only_sinchoのみ適用可能。CNSを今後適用予定")
            st.write("ChemTSの設定が完了しました。上記の設定で問題無ければ、次のステップに進んでください。")
            st.success("ChemTSの設定が完了しました。上記の設定で問題無ければ、次のステップに進んでください。")

        if sub_tab == "AAScore Settings":
            st.title("AAScore Settings")
            if "aascore_settings" not in st.session_state:
                st.session_state.aascore_settings = {}
            if "method" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["method"] = "all"  # スコア計算の方法
            if "num_of_cpd" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["num_of_cpd"] = 50  # ランダム選択する化合物数
            if "reward_cutoff" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["reward_cutoff"] = 1.00  # カットオフのreward値
            if "conf_per_cpd" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["conf_per_cpd"] = 20  # 1化合物当たりのconformation数
            if "max_attempts" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["max_attempts"] = 100  # Embedの最大試行回数
            if "rms_thresh" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["rms_thresh"] = 0.25  # 構造間でpruneするRMS閾値[Å]
            if "protein_range" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["protein_range"] = 13  # 計算時のポケット範囲(ヒットからX[Å]以内の残基)
            if "output_num" not in st.session_state.aascore_settings:
                st.session_state.aascore_settings["output_num"] = 5000   # 出力sdfファイルに格納する化合物数



            st.session_state.aascore_settings["method"] = st.selectbox("生成された化合物の内、どれをスコア計算するか？", options=["all", "rand"], index=["all", "rand"].index(st.session_state.aascore_settings["method"]))
            if st.session_state.aascore_settings["method"] == "rand":
                st.session_state.aascore_settings["num_of_cpd"] = st.number_input("ランダムに選択する化合物の数", value=st.session_state.aascore_settings["num_of_cpd"], step=1)
            st.session_state.aascore_settings["reward_cutoff"] = st.number_input("カットオフのreward値", value=st.session_state.aascore_settings["reward_cutoff"], step=0.01, max_value=1.00, min_value=0.00)
            st.session_state.aascore_settings["conf_per_cpd"] = st.number_input("1化合物当たりの最大conformation数", value= st.session_state.aascore_settings["conf_per_cpd"], step=1)
            with st.expander("conformation生成の追加設定"):
                st.session_state.aascore_settings["max_attempts"] = st.number_input("Embedの最大試行回数", value= st.session_state.aascore_settings["max_attempts"], step=1)
                st.session_state.aascore_settings["rms_thresh"] = st.number_input("構造間でpruneするRMS閾値[Å]", value= st.session_state.aascore_settings["rms_thresh"], step=0.01)
            st.session_state.aascore_settings["protein_range"] = st.number_input("計算時のポケット範囲(ヒットからX[Å]以内の残基)", value= st.session_state.aascore_settings["protein_range"], step=1)
            st.session_state.aascore_settings["output_num"] = st.number_input("出力sdfファイルに格納する化合物数", value= st.session_state.aascore_settings["output_num"], step=1)
            st.success("AAScoreの設定が完了しました。上記の設定で問題無ければ、次のステップに進んでください。")

        if sub_tab == "Summary":

            param_lines = ""
            for addparam in st.session_state.md_settings["additional_parameters"]:
                param_lines += "- "+os.path.join(st.session_state.general_settings["directory"], "99_TMP", addparam)+"\n      "

            replace_dict = {
                        "__OUTDIR__": str(st.session_state.general_settings["directory"]),
                        "__NUM_THREADS__": str(st.session_state.general_settings["use_num_threads"]),
                        "__INPUT_COMPLEX__": str(os.path.join(st.session_state.general_settings["directory"], "99_TMP", st.session_state.uploaded_pdb_file.name)),
                        "__HIT_RESNAME__": str(st.session_state.hit_residue.split(" ")[0]),
                        "__FORCE_FIELD_PROTEIN__": str(st.session_state.md_settings["force_field"][0]),
                        "__FORCE_FIELD_LIGAND__": str(st.session_state.md_settings["force_field"][1]),
                        "__FORCE_FIELD_WATER__": str(st.session_state.md_settings["force_field"][2].lower()),
                        "#__ADDITIONAL_PARAMS__": param_lines,
                        "__BOX_SHAPE__": str(st.session_state.md_settings["box_shape"]),
                        "__BOX_SIZE__": str(st.session_state.md_settings["box_size"]),
                        "__BUFFER__": str(st.session_state.md_settings["buffer"]),
                        # "__TEMPERATURE__": str(st.session_state.md_settings["temperature"]),
                        "__PRODUCTION_RUNTIME__": str(int(st.session_state.md_settings["pr_run_time"])*1000),
                        "__PRODUCTION_REC_INTERVAL__": str(int(st.session_state.md_settings["pr_rec_interval"])*1000),
                        "__SNAPSHOTS__": str(st.session_state.md_settings["snapshots"]),
                        "__SINCHO_DISTANCE_RANGE__": str(st.session_state.p2c_sincho_settings["distance_range"]),
                        "__SINCHO_NPAIRS_PER_SNAP__": str(st.session_state.p2c_sincho_settings["npairs_per_snap"]),
                        "__SINCHO_FOR_CHEMTS__": str(st.session_state.p2c_sincho_settings["for_chemts"]),
                        "__SINCHO_RESTRICT_RPOINTS__": str(st.session_state.p2c_sincho_settings["r_point_atoms"]),
                        "__CHEMTS_NUM_LOOPS__": str(st.session_state.chemts_settings["num_chemts_loops"]),
                        "__CHEMTS_C_VAL__": str(st.session_state.chemts_settings["c_val"]),
                        "__CHEMTS_THRESHOLD_TYPE__": str(st.session_state.chemts_settings["threshold_type"]),
                        "__CHEMTS_THRESHOLD__": str(st.session_state.chemts_settings["threshold"]),
                        #"__CHEMTS_FUNCTION_FORMAT__": str(st.session_state.chemts_settings["function_format"]),
                        "__AASCORE_METHOD__": str(st.session_state.aascore_settings["method"]),
                        "__AASCORE_NUM_OF_CPD__": str(st.session_state.aascore_settings["num_of_cpd"]),
                        "__AASCORE_REWARD_CUTOFF__": str(st.session_state.aascore_settings["reward_cutoff"]),
                        "__AASCORE_CONF_PER_CPD__": str(st.session_state.aascore_settings["conf_per_cpd"]),
                        "__AASCORE_MAX_ATTEMPTS__": str(st.session_state.aascore_settings["max_attempts"]),
                        "__AASCORE_RMS_THRESH__": str(st.session_state.aascore_settings["rms_thresh"]),
                        "__AASCORE_PROTEIN_RANGE__": str(st.session_state.aascore_settings["protein_range"]),
                        "__AASCORE_OUTPUT_NUM__": str(st.session_state.aascore_settings["output_num"]),

                    }


            st.title("Settings Summary")
            st.write(f"ディレクトリ：{os.getcwd()}に入力ファイルを作成します。")
            st.session_state.yaml_name = st.text_input("YAMLファイル名", value="conditions_lala.yaml")
            if st.session_state.yaml_content is None:
                with open(os.path.join(os.path.dirname(__file__), "conditions_tmp.yaml"), 'r') as file:
                    yaml_content = file.read()
                    for k, v in replace_dict.items():
                        yaml_content = yaml_content.replace(k, v)
                st.session_state.yaml_content = yaml_content
                for k in replace_dict.keys():
                    if k in st.session_state.yaml_content:
                        st.warning(f"{k} が 指定されていません。他のタブで設定し直してください。")


            st.session_state.yaml_content = st.text_area(
                    "Edit YAML content: 追加で編集したい場合は以下を変更してください。",
                    value=st.session_state.yaml_content,
                    height=600,
                    key="edited_yaml"
                )
            if st.button("yamlの初期化"):
                st.session_state.yaml_content = None
                st.success("YAML内容が初期化されました。")
                try:
                    st.rerun()
                except Exception as e:
                    st.experimental_rerun()
            if st.button("Save YAML to File"):
                overwrite_avoider = 0
                yaml_file = os.path.join(os.getcwd(), st.session_state.yaml_name)
                if not yaml_file.endswith(".yaml"):
                    yaml_file += ".yaml"
                    st.session_state.yaml_name += ".yaml"
                # ファイル名が既に存在する場合は、連番を付けて保存
                while os.path.exists(yaml_file):
                    overwrite_avoider += 1
                    yaml_file = os.path.join(os.getcwd(), f"{st.session_state.yaml_name[:-5]}_{overwrite_avoider}.yaml")


                edited_text = st.session_state.yaml_content
                with open(yaml_file, 'w') as f:
                    f.write(edited_text)
                st.success(f"YAMLファイルを保存しました: {yaml_file}")
                if overwrite_avoider > 0:
                    st.warning(f"ファイル名が重複したため、ファイル名に'_{overwrite_avoider}'を付加しました。")
                try:
                    shutil.copy(yaml_file, os.path.join(st.session_state.general_settings['directory'], st.session_state.general_settings['tmp_dir']))
                    st.success(f"{yaml_file}を{os.path.join(st.session_state.general_settings['directory'], st.session_state.general_settings['tmp_dir'])}にバックアップしました。")
                except:
                    st.error(f"{yaml_file}の{os.path.join(st.session_state.general_settings['directory'], st.session_state.general_settings['tmp_dir'])}に対するバックアップに失敗しました。手動でバックアップしてください。")
                #extend_driver.pyがos.getcwd()になければ、コピーしてくる。
                if not os.path.exists(os.path.join(os.getcwd(), "extend_driver.py")):
                    try:
                        shutil.copy(os.path.join(os.path.dirname(__file__), "extend_driver.py"), os.getcwd())
                        st.success("extend_driver.pyを現在の作業ディレクトリにコピーしました。")
                        st.success("これで全ての設定が完了しました。以下のコマンドを実行して処理を開始してください。")
                        st.code(f"python extend_driver.py {yaml_file.split('/')[-1]}")
                    except Exception as e:
                        st.error(f"extend_driver.pyのコピーに失敗しました: {e}")
                else:
                    st.success("extend_driver.pyは既に現在の作業ディレクトリに存在します。")
                    st.success("これで全ての設定が完了しました。以下のコマンドを実行して処理を開始してください。")
                    st.code(f"python extend_driver.py {yaml_file.split('/')[-1]}")



    def _pdb_3dview(self, pdbfile, zoomres = None):
        #PDBファイルを3次元表示
        pdb_str = st.session_state.uploaded_pdb_file.getvalue().decode("utf-8")
        view = py3Dmol.view(height=500, width=800)
        view.addModel(pdb_str, "pdb")

        # 水分子（HOH, WAT）
        view.setStyle({'resn': 'HOH'}, {"sphere": {"color": "skyblue", "radius": 0.3}})
        view.setStyle({'resn': 'WAT'}, {"sphere": {"color": "skyblue", "radius": 0.3}})
        # リガンドなど（水以外の HETATM）
        view.setStyle({'hetflag': True, 'resn': ['HOH', 'WAT'], 'invert': True},
                      {"stick": {"colorscheme": "greenCarbon"}})
        # タンパク質
        view.setStyle({'hetflag': False}, {"cartoon": {"color": "gray"}})
        if zoomres:
            view.zoomTo({'resi': zoomres.split(" ")[1]})
            view.setStyle({'resi':zoomres.split(" ")[1]}, {"stick":{"colorscheme": "orangeCarbon"}})
        else:
            view.zoomTo()
        html(view._make_html(), height=500, width=800)
        return
    
    def _pdb_3dview_res(self, res):
        #PDBファイルを3次元表示
        pdb_str = st.session_state.uploaded_pdb_file.getvalue().decode("utf-8")
        view = py3Dmol.view(height=500, width=800)
        view.addModel(pdb_str, "pdb")

        name_list = []
        view.setStyle({'resi': res.split(" ")[1]}, {"stick":{"colorscheme": "orangeCarbon"}})
        view.zoomTo({'resi': res.split(" ")[1]})
        for line in pdb_str.splitlines():
            if (line.startswith("ATOM") or line.startswith("HETATM")) and line[22:26].strip() == res.split(" ")[1] and "H" not in line[12:16].strip():
                name = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                view.addLabel(name, {
                    "position": {"x": x, "y": y, "z": z},
                    "fontColor": "black",
                    "backgroundColor": "white",
                    "fontSize": 15,
                    "inFront": True
                })
                name_list.append(name)

        html(view._make_html(), height=500, width=800)
        return name_list
    
    def _pdb_3dview_multires(self, resseq=list):
        #PDBファイルを3次元表示
        pdb_str = st.session_state.uploaded_pdb_file.getvalue().decode("utf-8")
        view = py3Dmol.view(height=500, width=800)
        view.addModel(pdb_str, "pdb", {"keepH": True})
        view.setStyle({'hetflag': False}, {"cartoon": {"color": "gray"}})
        view.zoomTo()

        for r in resseq:
            xyz = []
            xyz = [[float(line[30:38]),float(line[38:46]),float(line[46:54])] for line in pdb_str.splitlines() if (line[22:26].strip() == str(r) and (line.startswith("ATOM  ") or line.startswith("HETATM")))] 
            x = sum([i[0] for i in xyz])/len(xyz)
            y = sum([i[1] for i in xyz])/len(xyz)
            z = sum([i[2] for i in xyz])/len(xyz)
            view.setStyle({'resi': r}, {"stick":{"colorscheme": "greenCarbon"}, "hydrogens": True})
            view.addLabel(r, {
                    "position": {"x": x, "y": y, "z": z},
                    "fontColor": "black",
                    "backgroundColor": "white",
                    "fontSize": 15,
                    "inFront": True
                })
        wrapped_html = f"""
            <div style="width:100%; height:100%; overflow:hidden;">
                {view._make_html()}
            </div>
            """
        #html(view._make_html(), height=500, scrolling=True)
        html(wrapped_html, height=500, scrolling=True)
        return 
    

    def _residue_parser(self, pdbfile):
        # PDB ファイルを解析して残基一覧を取得
        parser = PDBParser(QUIET=True)
        # BytesIO オブジェクトをテキスト形式に変換
        pdb_text = io.StringIO(pdbfile.getvalue().decode("utf-8"))
        structure = parser.get_structure("uploaded_structure", pdb_text)
        reslist = [
            f"{residue.get_resname()} {residue.id[1]}"
            for model in structure
            for chain in model
            for residue in chain
        ]
        return reslist


    def _pdb_editable_board(self):
        st.markdown("### 🧬 PDB Structure Refinement / Editing Panel")
        # ============================================
        # Load and parse PDB text
        # ============================================

        #
        if st.session_state.edited_pdb_content is None:
            pdb_str = st.session_state.uploaded_pdb_file.getvalue().decode("utf-8")
        else:
            pdb_str = st.session_state.edited_pdb_content.getvalue().decode("utf-8")
        ligand_resname = st.session_state.hit_residue.split(" ")[0]  # ex) "LIG"

        df, other_lines = self._parse_pdb(pdb_str)
        st.divider()

        replaces_dict_cys, replaces_dict_his, replaces_dict_liganame, replaces_dict_protonated = {}, {}, {}, {}

        # ============================================
        # 1. Cysteine State Handling
        with st.expander("Cysteine definition (CYS / CYX / CYM)", expanded=True):
            cys_like = df[df["resname"].str.strip().isin(["CYS", "CYX", "CYM"])]

            if cys_like.empty:
                st.info("No cysteine-type residues (CYS/CYX/CYM) found.")
            else:
                # 残基リスト（重複なし）
                cys_res = (
                    cys_like[["chain", "resname", "resseq", "icode"]]
                    .drop_duplicates()
                    .sort_values(["chain", "resseq", "icode"]))

                # 3D可視化：対象残基をハイライト表示
                # （ここは既存の関数に残基番号一覧を渡すだけでOK）
                self._pdb_3dview_multires(cys_res["resseq"].tolist())

                st.markdown("#### Cysteine state per residue")
                st.write('you specify the cysteine states')
                st.write('SH State->"CYS", S-S State->"CYX", S- State->"CYM"')

                # GUI で CYS/CYX/CYM を選択させる
                state_map = {}
                for _, r in cys_res.iterrows():
                    chain = r["chain"]
                    resseq = r["resseq"]
                    icode = r["icode"]
                    current = r["resname"].strip()  # CYS/CYX/CYM
                    resid_label = f"Chain {chain}, Residue {resseq}{icode.strip() or ''} ({current})"
                    options = ["CYS", "CYX", "CYM"]
                    default_idx = options.index(current) if current in options else 0

                    selected = st.selectbox(
                        resid_label,
                        options=options,
                        index=default_idx,
                        key=f"cys_state_{chain}_{resseq}_{icode}")
                    # キー: (chain, resseq, icode) → 値: 選択されたstate
                    state_map[(current, chain, resseq, icode)] = selected

                for (orig, chain, resseq, icode), new_state in state_map.items():
                    if orig != new_state:
                        replaces_dict_cys[f"{orig} {chain}{str(resseq).rjust(4)}"] = f"{new_state} {chain}{str(resseq).rjust(4)}"
        st.divider()
        # ============================================
        # 2. Histidine Protonation (HID/HIE/HIP)
        """
        with st.expander("2. Histidine protonation (HID / HIE / HIP)", expanded=False):
            his_like = df[df["resname"].str.strip().isin(["HIS", "HID", "HIE", "HIP"])]

            if his_like.empty:
                st.info("No histidine-type residues (HIS/HID/HIE/HIP) found.")
            else:
                his_res = (
                    his_like[["chain", "resname", "resseq", "icode"]]
                    .drop_duplicates()
                    .sort_values(["chain", "resseq", "icode"]))

                # 3D可視化：対象残基をハイライト表示
                # （ここは既存の関数に残基番号一覧を渡すだけでOK）
                self._pdb_3dview_multires(his_res["resseq"].tolist())

                st.markdown("#### Histidine state per residue")
                st.write('you specify the histidine states')
                st.write('δ State->"HID", ε State->"HIE", Protonated State->"HIP"')

                state_map = {}
                for _, r in his_res.iterrows():
                    chain = r["chain"]
                    resseq = r["resseq"]
                    icode = r["icode"]
                    current = r["resname"].strip()  # CYS/CYX/CYM
                    resid_label = f"Chain {chain}, Residue {resseq}{icode.strip() or ''} ({current})"
                    options = ["HIS", "HID", "HIE", "HIP"]
                    default_idx = options.index(current) if current in options else 0

                    selected = st.selectbox(
                        resid_label,
                        options=options,
                        index=default_idx,
                        key=f"his_state_{chain}_{resseq}_{icode}")
                    # キー: (chain, resseq, icode) → 値: 選択されたstate
                    state_map[(current, chain, resseq, icode)] = selected
                for (orig, chain, resseq, icode), new_state in state_map.items():
                    if orig != new_state:
                        replaces_dict_his[f"{orig} {chain}{str(resseq).rjust(4)}"] = f"{new_state} {chain}{str(resseq).rjust(4)}"
        st.divider()
        """

        # ============================================
        # other residues with protonated states
        with st.expander("Other residues with protonated states", expanded=False):
            st.info("ASP or ASH")
            asp_like = df[df["resname"].str.strip().isin(["ASP"])]
            #aspの残基情報(["chain"], ["resseq"])の一覧を取得
            chain_resseq_unique = asp_like.drop_duplicates(subset=["chain", "resseq"])
            for _, r in chain_resseq_unique.iterrows():
                chain = r["chain"]
                resseq = r["resseq"]
                #chainとresseqでフィルタリングして、その残基の原子名一覧を取得し、"HD2"があればst.write()で表示
                atom_names = asp_like[(asp_like["chain"] == chain) & (asp_like["resseq"] == resseq)]["name"].tolist()
                if "HD2" in [ an.strip() for an in atom_names]:
                    replaces_dict_protonated[f"ASP {chain}{str(resseq).rjust(4)}"] = f"ASH {chain}{str(resseq).rjust(4)}"
                    st.write(f"Residue ASP {chain}{resseq} is protonated to ASH.")
            st.info("GLU or GLH")
            glu_like = df[df["resname"].str.strip().isin(["GLU"])]
            #gluの残基情報(["chain"], ["resseq"])の一覧を取得
            chain_resseq_unique = glu_like.drop_duplicates(subset=["chain", "resseq"])
            for _, r in chain_resseq_unique.iterrows():
                chain = r["chain"]
                resseq = r["resseq"]
                #chainとresseqでフィルタリングして、その残基の原子名一覧を取得し、"HE2"があればst.write()で表示
                atom_names = glu_like[(glu_like["chain"] == chain) & (glu_like["resseq"] == resseq)]["name"].tolist()
                if "HE2" in [ an.strip() for an in atom_names]:
                    replaces_dict_protonated[f"GLU {chain}{str(resseq).rjust(4)}"] = f"GLH {chain}{str(resseq).rjust(4)}"
                    st.write(f"Residue GLU {chain}{resseq} is protonated to GLH.")

            st.info("HID or HIE or HIP")
            his_like = df[df["resname"].str.strip().isin(["HIS"])]
            #aspの残基情報(["chain"], ["resseq"])の一覧を取得
            chain_resseq_unique = his_like.drop_duplicates(subset=["chain", "resseq"])
            for _, r in chain_resseq_unique.iterrows():
                chain = r["chain"]
                resseq = r["resseq"]
                #chainとresseqでフィルタリングして、その残基の原子名一覧を取得し、"HD2"があればst.write()で表示
                atom_names = his_like[(his_like["chain"] == chain) & (his_like["resseq"] == resseq)]["name"].tolist()
                atom_names = [an.strip() for an in atom_names]
                if "HD2" in atom_names and "HE2" in atom_names:
                    replaces_dict_protonated[f"HIS {chain}{str(resseq).rjust(4)}"] = f"HIP {chain}{str(resseq).rjust(4)}"
                    st.write(f"Residue HIS {chain}{resseq} is protonated to HIP.")
                elif "HD2" in atom_names and "HE2" not in atom_names:
                    replaces_dict_protonated[f"HIS {chain}{str(resseq).rjust(4)}"] = f"HID {chain}{str(resseq).rjust(4)}"
                    st.write(f"Residue HIS {chain}{resseq} is protonated to HID.")
                elif "HD2" not in atom_names and "HE2" in atom_names:
                    replaces_dict_protonated[f"HIS {chain}{str(resseq).rjust(4)}"] = f"HIE {chain}{str(resseq).rjust(4)}"
                    st.write(f"Residue HIS {chain}{resseq} is protonated to HIE.")


        # ============================================
        # 3. ACE / NME capping residue validation
        with st.expander("5. Capping residues (ACE / NME)", expanded=False):
            replaces_dict_cap, delete_list = self._check_capping(df)
        st.divider()
        # ============================================
        # 4. Ligand atom name duplication check
        with st.expander("4. Ligand atom name duplication check", expanded=False):
            dup = self._check_dup_atom_names(df, ligand_resname)
            if dup.empty:
                st.success("✔ No duplicated atom names detected.")
            else:
                
                st.error("❗Atom name duplication detected in the ligand:")
                st.write("Current")
                st.dataframe(dup)
                st.warning("Refine the atom names in the table below (within 4 characters)")
                df = st.data_editor(dup, use_container_width=True, key="atom_editor")

                # replace dict
                for (idx_bef, bef), (idx_aft, aft) in zip(dup.iterrows(), df.iterrows()):
                    if bef["name"] != aft["name"]:
                        b = f'{bef["serial"]} {bef["name"].strip().rjust(4)}'
                        a = f'{aft["serial"]} {aft["name"].strip().rjust(4)}'
                        replaces_dict_liganame[b] = a
        st.divider()

        delete_list.append("CONECT")#CONECT行があると金属がBonded扱いになって煩雑になるので。



        # ============================================
        """
        # 5. Selection of A/B terminal residues (XXXA / XXXB) Not Yet Implemented
        with st.expander("3. A / B chain residue selection (terminal variants)", expanded=False):
            # TODO: 実装例: residue options
            st.info("Select A/B variants (coming soon).")
        st.divider()
        """
        # ============================================

        # ============================================
        # 6. Output final edited PDB
        # ============================================
        #replaces群をpdbファイルに置換させる
        bef_pdb = pdb_str.splitlines()
        for reps in [replaces_dict_cys, replaces_dict_his, replaces_dict_protonated, replaces_dict_cap, replaces_dict_liganame]:
            for rep_b, rep_a in reps.items():
                # rep_b: 置換前の文字列、rep_a: 置換後の文字列でsplitlinesごとに置換
                for i, line in enumerate(bef_pdb):
                    if rep_b in line:
                        bef_pdb[i] = line.replace(rep_b, rep_a)
        for i, line in enumerate(bef_pdb):
            for del_str in delete_list:
                if del_str in line:
                    bef_pdb[i] = ""
        bef_pdb = [line for line in bef_pdb if line != ""]  # 空行削除
        #新しいpdbファイルとして保存
        pdb_str_edited = "\n".join(bef_pdb)
        #このままpdbファイルとして保存する

        st.button("Apply Revisions and Generate Edited PDB", key="apply_pdb_edits")
        if st.session_state.get("apply_pdb_edits"):
            tmp_path = os.path.join(st.session_state.general_settings["tmp_dir"], "edited_input.pdb")
            with open(tmp_path, "w") as f:
                f.write(pdb_str_edited)
            st.success(f"Edited PDB file has been generated and saved to {tmp_path}.")
            edited_pdb_bytes = pdb_str_edited.encode("utf-8")
            buf = io.BytesIO(edited_pdb_bytes)
            buf.name = "edited_input.pdb"
            st.session_state.edited_pdb_content = buf
            st.session_state.uploaded_pdb_file = st.session_state.edited_pdb_content

        """
        st.download_button(
            "⬇ Download edited PDB",
            data=pdb_str_edited,
            file_name="edited2.pdb",
            mime="chemical/x-pdb"
        )
        """



    # ================================
    # 1. PDBパース & 再構築
    # ================================
    def _parse_pdb(self, pdb_text: str):
        """
        PDB文字列を Atom/HETATM の DataFrame と その他行のlist に分解
        """
        lines = pdb_text.splitlines()
        atom_records = []
        other_lines = []

        for line in lines:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 54:
                atom_records.append({
                    "record": line[0:6],               # "ATOM  " / "HETATM"
                    "serial": int(line[6:11]),
                    "name": line[12:16],               # そのまま保持（stripは表示側で）
                    "altloc": line[16],
                    "resname": line[17:20],
                    "chain": line[21],
                    "resseq": int(line[22:26]),
                    "icode": line[26],
                    "x": float(line[30:38]),
                    "y": float(line[38:46]),
                    "z": float(line[46:54]),
                    "occupancy": line[54:60],
                    "tempfactor": line[60:66],
                    "segment": line[72:76] if len(line) >= 76 else "    ",
                    "element": line[76:78] if len(line) >= 78 else "  ",
                    "charge": line[78:80] if len(line) >= 80 else "  ",
                })
            else:
                other_lines.append(line)

        df = pd.DataFrame(atom_records)
        return df, other_lines


    def _rebuild_pdb(self, df: pd.DataFrame, other_lines):
        """
        DataFrame + その他行 から PDB文字列を再構築
        """
        atom_lines = []
        # 並び順は適当に record→chain→resseq→serial
        for _, r in df.sort_values(["record", "chain", "resseq", "serial"]).iterrows():
            line = (
                f"{str(r['record']):6s}"
                f"{int(r['serial']):5d} "
                f"{str(r['name']):<4s}"
                f"{str(r['altloc']):1s}"
                f"{str(r['resname']):>3s} "
                f"{str(r['chain']):1s}"
                f"{int(r['resseq']):4d}"
                f"{str(r['icode']):1s}   "
                f"{float(r['x']):8.3f}"
                f"{float(r['y']):8.3f}"
                f"{float(r['z']):8.3f}"
                f"{str(r['occupancy']):>6s}"
                f"{str(r['tempfactor']):>6s}      "
                f"{str(r['element']):>2s}"
                f"{str(r['charge']):>2s}"
            )
            atom_lines.append(line)

        # ENDが既にother_linesに含まれていても、最後にも付けてしまう簡易実装
        all_lines = other_lines + atom_lines + ["END"]
        return "\n".join(all_lines)

    # ================================
    # 4. リガンドAtomName重複チェック
    # ================================
    def _check_dup_atom_names(self, df: pd.DataFrame, lig_resname: str):
        """
        指定リガンド(resname)について、同一残基内でAtomNameが重複しているものを返す。
        戻り値: 重複行だけのDataFrame（なければ空DataFrame）
        """
        lig_mask = df["resname"].str.strip() == lig_resname.strip()
        lig_df = df[lig_mask].copy()
        if lig_df.empty:
            # 空のDataFrameを返す
            return lig_df

        # 同一残基( chain, resseq, icode ) 内での重複をチェック
        dup_rows = []

        for (chain, resseq, icode), sub in lig_df.groupby(["chain", "resseq", "icode"]):
            # atom nameで重複
            duplicated = sub[sub.duplicated(subset=["name"], keep=False)]
            if not duplicated.empty:
                dup_rows.append(duplicated)

        if dup_rows:
            return pd.concat(dup_rows, axis=0)
        else:
            return lig_df.iloc[0:0]  # 空DataFrame


    # ================================
    # 5. ACE / NME capping 残基チェック
    # ================================
    ############################## 一旦Maestro形式の入力を想定してつくるよ ################################
    def _check_capping(self, df: pd.DataFrame):
        #最後に文字列置き換えをするための辞書
        replaces_dict = {}
        delete_list = []

        #NMAエラーチェッカー
        if "NMA" in df["resname"].values:
            st.write("NME residue is named as NMA in the PDB file. This needs to rename to 'NME' ")
            replaces_dict["NMA"]="NME"

        #ACE原子ラベルチェッカー
        #dfから"resname"=="ACE"の原子行を抜き出したい
        ace_atomnames = df[df["resname"].str.strip() == "ACE"]["name"].tolist()
        maestro_to_amber_ace = {"1H   ACE":"HH31 ACE", "2H   ACE":"HH32 ACE", "3H   ACE":"HH33 ACE"}
        for name in ace_atomnames:
            for k,v in maestro_to_amber_ace.items():
                if name.strip() in k:
                    replaces_dict[k]=v
                    st.write(f"ACE atom name '{name.strip()}' needs to be renamed to '{v}'")
                    break
        #ACEの次の残基のH1チェッカー
        ace_resids = df[df["resname"].str.strip() == "ACE"][["chain", "resseq", "icode"]].drop_duplicates()
        for _, r in ace_resids.iterrows():
            chain, resseq, icode = r["chain"], r["resseq"], r["icode"]
            next_res_mask = (
                (df["chain"] == chain) &
                (df["resseq"] == resseq + 1) &
                (df["icode"] == icode)
            )
            next_res_df = df[next_res_mask]
            if not next_res_df.empty:
                next_res_h1 = next_res_df[next_res_df["name"].str.strip() == "H1"]
                if not next_res_h1.empty:
                    resname = next_res_df.iloc[0]["resname"].strip()
                    st.write(f"The H1 atom in the residue following ACE (residue {resname} {resseq + 1}) needs to be deleted to avoid duplication.")
                    delete_list.append("H1  "+resname)

        #NME原子ラベルチェッカー("NME" or "NMA")
        nme_atomnames = df[df["resname"].str.strip().isin(["NME", "NMA"])]["name"].tolist()
        #minicondaの実装ではCH3、
        maestro_to_amber_nme = {" CA  NME":" CH3 NME", "1HA  NME":"HH31 NME", "2HA  NME":"HH32 NME", "3HA  NME":"HH33 NME", }
        for name in nme_atomnames:
            for k,v in maestro_to_amber_nme.items():
                if name.strip() in k:
                    replaces_dict[k]=v
                    st.write(f"NME atom name '{name.strip()}' needs to be renamed to '{v}'")
                    break
    
        return replaces_dict, delete_list


    # ================================
    # 汎用ユーティリティ
    # ================================
    def _residue_id(self, row):
        """chain:resname:resseq(icode) のようなID文字列を作る"""
        icode = (row.get("icode", " ") or " ").strip()
        icode_str = icode if icode else ""
        return f"{row['chain']}:{row['resname'].strip()}:{row['resseq']}{icode_str}"






    


