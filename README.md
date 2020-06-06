random-signal-process-experiment
================================

A backup of my homework.

仿真线性调频脉冲雷达的信号处理。设线性调频带宽为10MHz，时宽为200us，占空比10%
，雷达载频为10GHz，输入噪声为高斯白噪声。目标模拟分单目标和双目标两种情况，目
标回波输入信噪比可变（-40dB $\sim$ 10dB），目标速度可变（0 $\sim$ 1000m/s），
目标距离可变（0 $\sim$ 10000\m），相干处理时间不小于100ms。

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
latexmk -pvc main.tex
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
