using Plots
using DifferentialEquations

N = 11700
I0 = 270 # заболевшие особи
R0 = 49 # особи с иммунитетом
S0 = N - I0 - R0 # здоровые, но восприимчивые особи
<<<<<<< HEAD
alpha = 0.4 # коэффициент заболеваемости
beta = 0.1 # коэффициент выздоровления
=======
alpha = 0.4 # коэффициент заболеваемости????
beta = 0.1 # коэффициент выздоровления????
>>>>>>> b3fe2954b1cda54fcb6d1ddefeea5afc35f70508

#I0 > I*
function ode_fn(du, u, p, t)
    S, I, R = u
    du[1] = -alpha*u[1]
    du[2] = alpha*u[1] - beta*u[2]
    du[3] = beta*I
end

v0 = [S0, I0, R0]
tspan = (0.0, 120.0)
prob = ODEProblem(ode_fn, v0, tspan)
sol = solve(prob, dtmax=0.05)
S = [u[1] for u in sol.u]
I = [u[2] for u in sol.u]
R = [u[3] for u in sol.u]
T = [t for t in sol.t]

plt = plot(
  dpi=600,
  legend=:right)

plot!(
  plt,
  T,
  S,
  label="Восприимчивые особи",
  color=:blue)
plot!(
  plt,
  T,
  I,
  label="Инфицированные особи",
  color=:green)
plot!(
  plt,
  T,
  R,
  label="Особи с иммунитетом",
  color=:red)


savefig(plt, "lab06_2.png")