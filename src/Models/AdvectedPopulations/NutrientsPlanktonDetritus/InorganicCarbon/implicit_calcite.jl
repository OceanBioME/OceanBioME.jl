struct CarbonateSystem{N} <: AbstractInorganicCarbon end

function CarbonateSystem(replicates = 1)
    manifest_carbonate_replicates!(replicates)

    return CarbonateSystem{replicates}()
end

required_biogeochemical_tracers(::CarbonateSystem{1}) = (:DIC, :Alk)
required_biogeochemical_tracers(::CarbonateSystem{N}) where N = (map(n->Symbol(:DIC, n), 1:N)..., map(n->Symbol(:Alk, n), 1:N)...)
required_biogeochemical_auxiliary_fields(::CarbonateSystem) = tuple()

const NPD_DIC_Alk{FT} = NPD{FT, <:Any, <:Any, <:Any, <:CarbonateSystem}

function manifest_carbonate_replicates!(N)
    if N>1
      for n in 1:N
          DIC_name = Symbol(:DIC, n)
          Alk_name = Symbol(:Alk, n)
          @eval begin
              @inline (bgc::NPD_DIC_Alk)(i, j, k, grid, ::Val{$(QuoteNode(DIC_name))}, clock, fields, auxiliary_fields) =
                  bgc(i, j, k, grid, Val(:DIC), clock, fields, auxiliary_fields)
              @inline (bgc::NPD_DIC_Alk)(i, j, k, grid, ::Val{$(QuoteNode(Alk_name))}, clock, fields, auxiliary_fields) =
                  bgc(i, j, k, grid, Val(:Alk), clock, fields, auxiliary_fields)
          end
      end
    end
end

Base.summary(carbonates::CarbonateSystem{1}) = 
    string("CarbonateSystem $(required_biogeochemical_tracers(carbonates))")

Base.summary(carbonates::CarbonateSystem{N}) where N = 
    string("CarbonateSystem{realisations = $N} $(required_biogeochemical_tracers(carbonates))")

function Base.show(io::IO, c::CarbonateSystem{N}) where N
    msg = "CarbonateSystem $(required_biogeochemical_tracers(c))"

    if N>1
         msg *= "\n└── Realisations: $N"
    end

    print(io, msg)

    return nothing
end
