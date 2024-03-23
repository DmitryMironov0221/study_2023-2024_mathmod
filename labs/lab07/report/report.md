---
## Front matter
title: "Отчёт по лабораторной работе №7


Математическое моделирование"
subtitle: "Модель распространения рекламы. Вариант №32"
author: "Выполнил: Миронов Дмитрий Андреевич, 


НПИбд-02-21, 1032211701"



## Generic otions
lang: ru-RU
toc-title: "Содержание"

## Bibliography
bibliography: bib/cite.bib
csl: pandoc/csl/gost-r-7-0-5-2008-numeric.csl

## Pdf output format
toc: true # Table of contents
toc-depth: 2
lof: true # List of figures
fontsize: 12pt
linestretch: 1.5
papersize: a4
documentclass: scrreprt
## I18n polyglossia
polyglossia-lang:
  name: russian
  options:
	- spelling=modern
	- babelshorthands=true
polyglossia-otherlangs:
  name: english
## I18n babel
babel-lang: russian
babel-otherlangs: english
## Fonts
mainfont: PT Serif
romanfont: PT Serif
sansfont: PT Sans
monofont: PT Mono
mainfontoptions: Ligatures=TeX
romanfontoptions: Ligatures=TeX
sansfontoptions: Ligatures=TeX,Scale=MatchLowercase
monofontoptions: Scale=MatchLowercase,Scale=0.9
## Biblatex
biblatex: true
biblio-style: "gost-numeric"
biblatexoptions:
  - parentracker=true
  - backend=biber
  - hyperref=auto
  - language=auto
  - autolang=other*
  - citestyle=gost-numeric
## Pandoc-crossref LaTeX customization
figureTitle: "Рис."
tableTitle: "Таблица"
listingTitle: "Листинг"
lofTitle: "Список иллюстраций"
lolTitle: "Листинги"
## Misc options
indent: true
header-includes:
  - \usepackage{indentfirst}
  - \usepackage{float} # keep figures where there are in the text
  - \floatplacement{figure}{H} # keep figures where there are in the text
---

# Цель работы

Изучить и построить модель эффективности рекламы.

# Теоретическое введение. Построение математической модели.

Организуется рекламная кампания нового товара или услуги. Необходимо, чтобы прибыль будущих продаж с избытком покрывала издержки на рекламу. Вначале расходы могут превышать прибыль, поскольку лишь малая часть потенциальных покупателей будет информирована о новинке. Затем, при увеличении числа продаж, возрастает и прибыль, и, наконец, наступит момент, когда рынок насытиться, и рекламировать товар станет бесполезным.

Предположим, что торговыми учреждениями реализуется некоторая продукция, о которой в момент времени $t$ из числа потенциальных покупателей $N$ знает лишь $n$ покупателей. Для ускорения сбыта продукции запускается реклама по радио, телевидению и других средств массовой информации. После запуска рекламной кампании информация о продукции начнет распространяться среди потенциальных покупателей путем общения друг с другом. Таким образом, после запуска рекламных объявлений скорость изменения числа знающих о продукции людей пропорциональна как числу знающих о товаре покупателей, так и числу покупателей о нем не знающих

Модель рекламной кампании описывается следующими величинами.
Считаем, что $\frac{dn}{dt}$ - скорость изменения со временем числа потребителей, узнавших о товаре и готовых его купить,
$t$ - время, прошедшее с начала рекламной кампании,
$N$ - общее число потенциальных платежеспособных покупателей,
$n(t)$ - число  уже информированных клиентов.
Эта величина пропорциональна числу покупателей, еще не знающих о нем, это описывается следующим образом
$\alpha _1(t)(N-n(t))$, где $\alpha _1>0$ -  характеризует интенсивность рекламной кампании (зависит от затрат на рекламу в данный момент времени).
Помимо этого, узнавшие о товаре потребители также распространяют полученную информацию среди потенциальных покупателей, не знающих о нем (в этом случае работает т.н. сарафанное радио). Этот вклад в рекламу описывается величиной  $\alpha _2(t)n(t)(N-n(t))$. эта величина увеличивается с увеличением потребителей узнавших о товаре.

Математическая модель распространения рекламы описывается уравнением:

$$\frac{dn}{dt} = (\alpha _1(t) + \alpha _2(t)n(t))(N-n(t))$$

При $\alpha _1(t) >> \alpha _2(t)$ получается модель типа модели Мальтуса, решение которой имеет вид 

![График решения уравнения модели Мальтуса](image/00.PNG){ #fig:0001 width=70% height=70% }

В обратном случае $\alpha _1(t) << \alpha _2(t)$ получаем уравнение логистической кривой

![График логистической кривой](image/000.PNG){ #fig:0002 width=70% height=70% }

# Задание

**Вариант 32**

Постройте график распространения рекламы, математическая модель которой описывается следующим уравнением:

1.	$\frac{dn}{dt} = (0.66 + 0.000061n(t))(N-n(t))$
2.	$\frac{dn}{dt} = (0.000056 + 0.66n(t))(N-n(t))$
3.	$\frac{dn}{dt} = (0.66\sin{t} + 0.66\sin{(6t)}  n(t))(N-n(t))$

При этом объем аудитории $N = 860$, в начальный момент о товаре знает 2 человек.

Для случая 2 определите в какой момент времени скорость распространения рекламы будет иметь максимальное значение.

# Выполнение лабораторной работы

## Решение с помощью программ

### Julia

Код программы для первого случая $\frac{dn}{dt} = (0.66 + 0.000061n(t))(N-n(t))$:

```
using Plots
using DifferentialEquations

N = 860
n0 = 2

function ode_fn(du, u, p, t)
    (n) = u
    du[1] = (0.66 + 0.00061*u[1])*(N - u[1])
end

v0 = [n0]
tspan = (0.0, 30.0)
prob = ODEProblem(ode_fn, v0, tspan)
sol = solve(prob, dtmax = 0.05)
n = [u[1] for u in sol.u]
T = [t for t in sol.t]

plt = plot(
  dpi = 600,
  title = "Эффективность распространения рекламы (1) ",
  legend = false)
plot!(
  plt,
  T,
  n,
  color = :red)

savefig(plt, "lab07_1.png")
```

Код программы для второго случая $\frac{dn}{dt} = (0.000056 + 0.66n(t))(N-n(t))$:

```
using Plots
using DifferentialEquations

N = 860
n0 = 2

function ode_fn(du, u, p, t)
    (n) = u
    du[1] = (0.000056 + 0.66*u[1])*(N - u[1])
end

v0 = [n0]
tspan = (0.0, 0.1)
prob = ODEProblem(ode_fn, v0, tspan)
sol = solve(prob)
n = [u[1] for u in sol.u]
T = [t for t in sol.t]

max_dn = 0;
max_dn_t = 0;
max_dn_n = 0;
for (i, t) in enumerate(T)
    if sol(t, Val{1})[1] > max_dn
        global max_dn = sol(t, Val{1})[1]
        global max_dn_t = t
        global max_dn_n = n[i]
    end
end

plt = plot(
  dpi = 600,
  title = "Эффективность распространения рекламы (2) ",
  legend = false)
plot!(
  plt,
  T,
  n,
  color = :red)
plot!(
  plt,
  [max_dn_t],
  [max_dn_n],
  seriestype = :scatter,
  color = :red)

savefig(plt, "lab07_2.png")
```

Код программы для третьего случая $\frac{dn}{dt} = (0.66\sin{t} + 0.66\sin{(6t)}  n(t))(N-n(t))$:

```
using Plots
using DifferentialEquations

N = 860
n0 = 2

function ode_fn(du, u, p, t)
    (n) = u
    du[1] = (0.66 + 0.66*sin(6*t)*u[1])*(N - u[1])
end

v0 = [n0]
tspan = (0.0, 0.1)
prob = ODEProblem(ode_fn, v0, tspan)
sol = solve(prob, dtmax = 0.05)
n = [u[1] for u in sol.u]
T = [t for t in sol.t]

plt = plot(
  dpi = 600,
  title = "Эффективность распространения рекламы (3) ",
  legend = false)
plot!(
  plt,
  T,
  n,
  color = :red)

savefig(plt, "lab07_3.png")
```

### Результаты работы кода на Julia

![График распространения рекламы для первого случая, построенный на языке Julia](image/1.PNG){ #fig:001 width=70% height=70% }

![График распространения рекламы для второго случая, построенный на языке Julia](image/2.PNG){ #fig:002 width=70% height=70% }

![График распространения рекламы для третьего случая, построенный на языке Julia](image/3.PNG){ #fig:003 width=70% height=70% }

## OpenModelica

Код программы для первого случая $\frac{dn}{dt} = (0.66 + 0.000061n(t))(N-n(t))$:

```
model lab07_1
Real N = 610;
Real n;
initial equation
n = 2;
equation
der(n) = (0.66 + 0.00061*n)*(N-n);
end lab07_1;
```

Код программы для второго случая $\frac{dn}{dt} = (0.000056 + 0.66n(t))(N-n(t))$:

```
model lab07_2
Real N = 860;
Real n;
initial equation
n = 2;
equation
der(n) = (0.000056 + 0.66*n)*(N-n);
end lab07_2;
```

Код программы для третьего случая $\frac{dn}{dt} = (0.66\sin{t} + 0.66\sin{(6t)}  n(t))(N-n(t))$:

```
model lab07_3
Real N = 860;
Real n;
initial equation
n = 2;
equation
der(n) = (0.66 + 0.66*sin(6*time)*n)*(N-n);
end lab07_3;
```

### Результаты работы кода на OpenModelica

![График распространения рекламы для первого случая, построенный с помощью OpenModelica](image/4.PNG){ #fig:004 width=70% height=70% }

![График распространения рекламы для второго случая, построенный с помощью OpenModelica](image/5.PNG){ #fig:005 width=70% height=70% }

![График распространения рекламы для третьего случая, построенный с помощью OpenModelica](image/6.PNG){ #fig:006 width=70% height=70% }

# Анализ полученных результатов. Сравнение языков.

В итоге проделанной работы мы построили графики распространения рекламы для трех случаев на языках Julia и OpenModelica. Построение модели распространения рекламы на языке OpenModelica занимает значительно меньше строк, чем аналогичное построение на Julia. 
Кроме того, построения на языке OpenModelica проводятся относительно значения времени t по умолчанию, что упрощает нашу работу. 

# Вывод

В ходе выполнения лабораторной работы была изучена модель эффективности рекламы и в дальнейшем построена модель на языках Julia и Open Modelica.

# Список литературы. Библиография.

[1] Документация по Julia: https://docs.julialang.org/en/v1/

[2] Документация по OpenModelica: https://openmodelica.org/

[3] Решение дифференциальных уравнений: https://www.wolframalpha.com/

[4] Мальтузианская модель роста: https://www.stolaf.edu//people/mckelvey/envision.dir/malthus.html