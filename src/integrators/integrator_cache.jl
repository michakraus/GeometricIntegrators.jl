
abstract type IntegratorCache{DT,D} end

abstract type ODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type DAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type IDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PODEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end
abstract type PDAEIntegratorCache{DT,D} <: IntegratorCache{DT,D} end


IntegratorCache(params::Parameters) = error("IntegratorCache(params) not implemented for ", typeof(params))

create_integrator_cache(integrator::Integrator) = error("create_integrator_cache()! not implemented for ", typeof(integrator))
