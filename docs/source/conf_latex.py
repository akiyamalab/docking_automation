# LaTeX関連の設定を更新 (XeTeX用)
latex_engine = 'xelatex'
latex_elements = {
    'papersize': 'a4paper',
    'pointsize': '10pt',
    'preamble': r'''
\usepackage{xeCJK}
\setCJKmainfont{IPAMincho}
\setCJKsansfont{IPAGothic}
\setCJKmonofont{IPAGothic}
\XeTeXlinebreaklocale "ja"
\XeTeXlinebreakskip = 0pt plus 1pt
''',
    'figure_align': 'htbp',
    'fncychap': r'\usepackage[Bjarne]{fncychap}',
    'extraclassoptions': ',openany,oneside',
    'maketitle': r'\newcommand\sphinxbackoftitlepage{}\sphinxmaketitle',
    'sphinxsetup': r'warningBgColor={rgb}{1,1,1}',
}
