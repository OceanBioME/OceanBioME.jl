# [Sugar kelp (Saccharina latissima) individuals](@id SLatissima)

We have implemented a model of sugar kelp growth within this spatially infinitesimal Lagrangian particles framework originally based on the model of [Broch2012](@citet) and updated by [Broch2013](@citet), [Fossberg2018](@citet), and [Broch2019](@citet). This is the same model passively forced by [StrongWright2022](@citet).

The model tracks three variables, the frond area, A (dm²), carbon reserve, C (gC / gSW), and nitrate reserve, N (gN / gSW). The growth depends on the nitrate (and optionally ammonia) availability in the water, the temperature, and light availability. The minimum required coupling is with nitrates so the model can be coupled with an NPZD model, but can optionally uptake ammonia, DIC (CO₂), oxygen, and release dissolved organic matter (from exudation) and large detritus.  

Results could look something like this (from [StrongWright2022](@citet)):
![Example A, N, and C profiles from [StrongWright2022](@citet)](https://www.frontiersin.org/files/Articles/793977/fmars-08-793977-HTML/image_m/fmars-08-793977-g002.jpg)

## Model equations

As per [Broch2012](@citet) this model variables evolve as:

$\frac{dA}{dt} = \left(\mu - \nu\right)A,$

$\frac{dN}{dt} = J - \mu(N + N_\text{struct}),$

$\frac{dC}{dt} = P(1 - E) - R - \mu(C + C_\text{struct}).$

The apical frond loss given by:

$\nu = \frac{10^{-6}\exp(\varepsilon A)}{1 + 10^{-6}(\exp(\varepsilon A) - 1)}.$

The growth given by:

$\mu = f_\text{area}f_\text{seasonal}f_\text{temp}\min\left(\mu_c, \max(\mu_{NO_3}, \mu_{NH_4})\right),$

where:

$f_\text{area} = m_1\exp(-(A/A_0)^2) + m_2,$

$f_\text{seasonal} = a_1(1 + \text{sgn}(\lambda(n))|\lambda(n)|^{1/2}) + a_2,$

$f_\text{temp} = \left[ \begin{array}{ll}
                    0.08T + 0.2 & -1.8\leq T < 10 \\
                    1           & 10 \leq T \leq 15 \\
                    19/4 - T/4  & 15 < T \leq 19 \\
                    0           & T > 19
                 \end{array} \right]$

where $n$ is the day of the year, $\lambda$ is the normalised day length change, and $T$ is the temperature in degrees centigrade. The limiting rates ($\mu_c$, $\mu_{NO_3}$, $\mu_{NH_4}$) depend on the availability of carbon giving:

$\mu_c = 1 - \frac{C_\text{min}}{C},$

and on the available nitrogen which is either limited by the instantanious uptake of ammonia, or the nitrogen reserve. To find these limits $J$, the nutrient uptake, must first be found [Fossberg2018](@citep). The uptake is calculated by first finding the $NO_3$ uptake rate:

$J_{NO_3} = J_{NO_3\text{, max}}f_\text{curr}\frac{N_\text{max} - N}{N_\text{max} - N_\text{min}}\frac{NO_3}{k_{NO_3} + NO_3},$

where $f_{curr} = c_1(1 - \exp(-u / c_2)) + c_3$ [Broch2019](@citep) where $u$ is the relative current speed. We then calculate the potential ammonia uptake:

$\tilde{J}_{NH_4} = J_{NH_4\text{, max}}f_\text{curr}\frac{NH_4}{k_{NH_4} + NH_4}.$

This results in a theoretical instantaneous area increase rate of:

$\mu_{NH_4} = \frac{1}{N_\text{min} + N_\text{struct}}\frac{\tilde{J}_{NH_4}}{k_A},$

or growth from reserves of:

$\mu_{NO_3} = 1 - \frac{N_\text{min}}{N}.$

The actual resulting ammonia uptake is then:

$J_{NH_4} = \min\left(\tilde{J}_{NH_4}, \mu k_A (N_\text{min} + N_\text{struct})\right),$

and the total uptake $J = \frac{J_{NH_4} + J_{NO_3}}{k_A}$.

The production of carbon from photosynthesis is given by:

$P = P_S\left[1 - \exp\left(\frac{\alpha PAR}{P_S}\right)\right]\exp\left(\frac{\beta PAR}{P_S}\right),$

where $PAR$ is the photosynthetically available radiation and

$P_S = \frac{\alpha PAR}{\ln(1 + \alpha/\beta)}.$

$\alpha$ is a constant but $\beta$ depends on the maximum photosynthetic rate which is defined by both:

$P_\text{max} = \frac{\alpha PAR}{\ln(1 + \alpha/\beta)}\left(\frac{\alpha}{\alpha + \beta}\right)\left(\frac{\beta}{\alpha + \beta}\right) ^ {\beta/\alpha},$

and

$P_\text{max} = \frac{P_1\exp\left(\frac{T_{AP}}{T_{P1}} - \frac{T_{AP}}{T}\right)}{1 + \exp\left(\frac{T_{APL}}{T} - \frac{T_{APL}}{PL}\right) + \exp\left(\frac{T_{APH}}{T} - \frac{T_{APH}}{PH}\right)},$

where $T$ is the temperature in kelvin and the $T_X$ are Arrhenius temperature constants. We solve these iteratively to find $\beta$.

The exudation fraction is given by:

$E = 1 - \exp\left(\gamma(C_\min - C)\right).$

As per [Broch2013](@citet) $R$, the respiration rate, is given by:

$R = \left[R_A\left(\frac{\mu}{\mu_\text{max}} + \frac{J}{J_\text{max}}\right) + R_B\right]\exp\left(\frac{T_{ARR}}{T_1} - \frac{T_ARR}{T}\right)$

All other undefined symbols are constants which can be found in [Parameters](@ref parameters)

## Model conservations