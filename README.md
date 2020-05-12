random-signal-process-experiment
================================

A backup of my homework.

Dependent
---------

1.  A LaTex distribution. Such as [texlive].
2.  [pygments], syntax highlight.
3.  [boxie], syntax highlight.
4.  [njustthesis].
5.  [julia].

Install
-------

``` {.zsh}
git clone git@github.com:Freed-Wu/random-signal-process-experiment.git
cd random-signal-process-experiment
julia lst/main.jl # generate figures, tables, reference.
latexmk -pvc random-signal-process-experiment
```

Thanks
------

See `Thanks` in the paper.

Q & A
-----

More question see [Issues].

If you don't wanna complie, you can download the complied paper from
[Release]

  [texlive]: https://github.com/TeX-Live/texlive-source
  [pygments]: https://github.com/pygments/pygments
  [boxie]: https://github.com/registor/boxiesty
  [njustthesis]: https://github.com/Freed-Wu/njustthesis
  [julia]: https://github.com/JuliaLang/julia
  [Issues]: https://github.com/Freed-Wu/random-signal-process-experiment/issues
  [Release]: https://github.com/Freed-Wu/random-signal-process-experiment/releases/
