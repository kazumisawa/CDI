from cyvcf2 import VCF   # bgzip された VCF (.vcf.gz) を高速に読むためのライブラリ
import sys              # コマンドライン引数を受け取るため

# -----------------------------
# 入力の取得
# -----------------------------

vcf = VCF(sys.argv[1])  # 第1引数：VCFファイル（.vcf.gz）
region = sys.argv[2]    # 第2引数：解析する染色体領域（例 "chr1:100000-200000"）

# -----------------------------
# サンプル情報の取得
# -----------------------------

samples = vcf.samples   # VCF に含まれるサンプル ID のリスト
n = len(samples)        # サンプル数（個体数）

# -----------------------------
# CDI 計算に必要なカウンタ
# -----------------------------
# CDI = #(1/1) − 1/2 × #(0/1)
# なので、それぞれを数える

hom_alt = [0] * n       # ALT/ALT（1/1）：reference と異なるホモ接合サイト数
het     = [0] * n       # REF/ALT（0/1）：ヘテロ接合サイト数

# -----------------------------
# 指定領域内の variant を走査
# -----------------------------

for rec in vcf(region):         # 指定した染色体領域の各 variant を順に処理
    if len(rec.ALT) != 1:       # multi-allelic site（ALT が複数ある）
        continue                # CDI の理論は biallelic 前提なので無視

    # rec.genotypes は各サンプルについて
    # (a1, a2, phased) のタプルを返す
    #   a1, a2 : 0=REF, 1=ALT, -1=missing
    #   phased : True/False（今回は不要）

    for i, (a1, a2, _) in enumerate(rec.genotypes):
        if a1 < 0 or a2 < 0:    # missing genotype（./. など）
            continue            # 情報がないのでカウントしない

        if a1 == a2 == 1:       # 1/1（ALT ホモ接合）
            hom_alt[i] += 1     # #(1/1) をインクリメント

        elif a1 != a2:          # 0/1 または 1/0（ヘテロ接合）
            het[i] += 1         # #(0/1) をインクリメント

        # 0/0（REF ホモ接合）の場合は
        # CDI への寄与が 0 なので何もしない

# -----------------------------
# CDI の計算と出力
# -----------------------------

for i, s in enumerate(samples): # 各サンプル（個体）について
    # CDI = #(1/1) − 1/2 × #(0/1)
    cdi = hom_alt[i] - 0.5 * het[i]

    # サンプルIDとCDIを tab 区切りで出力
    print(f"{s}\t{cdi}")
