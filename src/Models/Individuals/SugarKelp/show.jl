import Base: summary, show

summary(::SugarKelp{FT}) where FT = "SugarKelp{FT} biogeochemistry"

show(io::IO, kelp::SugarKelp) = 
    print(io, summary(kelp), 
              " (Broch & Slagstad, 2012) tracking the `N`itrogen and `C`arbon in a frond of `A`rea")

