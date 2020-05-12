random-signal-process-experiment
================================

A backup of my homework.

Dependent
---------

1. A LaTex distribution. Such as [texlive].
2. [pygments](https://github.com/pygments/pygments), syntax highlight.
3. [boxie](https://github.com/registor/boxiesty), syntax highlight.
4. [njustthesis](https://github.com/Freed-Wu/njustthesis).
5. [julia](https://github.com/JuliaLang/julia).

Install
-------

``` {.zsh}
git clone git@github.com:Freed-Wu/random-signal-process-experiment.git
cd random-signal-process-experiment
julia lst/main.jl # generate figures, tables, reference.
latexmk -pvc random-signal-process-experiment
```

  [texlive]: https://github.com/TeX-Live/texlive-source
