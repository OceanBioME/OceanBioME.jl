struct CallableVal{V} end

CallableVal(x) = CallableVal{x}()

(::CallableVal{V})(args...; kwargs...) where {V} = V(args...; kwargs...)
