from setuptools import setup, find_packages


setup(
    name='humananalysis',
    version="0.0.12",                 # X.Y.Z 形式
    description="短めの説明",
    long_description="長めの説明",
    url='https://github.com/kyosukemine/humananalysis',
    author='kyosukemine',
    # author_email='メールアドレス',
    package_dir={'': 'humananalysis'},
    packages=find_packages('humananalysis', exclude=['usage_example']),
    install_requires=['numpy',
                      'pathlib2;python_version<"3.4"',
                      'scipy',
                      'gatspy',
                      'pywavelets',
                      'matplotlib',
                      'pandas'],
    license='GPL 3.0',
    classifiers=[
        # パッケージのカテゴリー
        # https://pypi.python.org/pypi?:action=list_classifiers
        # から該当しそなもを選んでおけばよい。
    ],
    keywords='キーワード',
    # install_requires=["依存関係のあるパッケージ"],
)
