export driftvelocity_point, driftvelocity_direction

"""
    driftvelocity_point(adf, target; normalize=false)
Evaluate the drift velocity of microbes towards a point `target`,
extracting their positions and velocities from the agent dataframe `adf`
(requires the existence of `:pos` and `:vel` fields).
Returns a matrix of instantaneous drift velocities with size
(nmicrobes, nsteps). By convention the drift velocity will be positive
for motion towards the target point.

If `normalize` is set to `true`, drift velocities are normalized by the
instantaneous speed of the microbe.
"""
function driftvelocity_point(adf, target; normalize=false)
    traj = vectorize_adf_measurement(adf, :pos)
    vels = vectorize_adf_measurement(adf, :vel)
    v_drift = zeros(size(traj)...)
    for k in eachindex(traj)
        pos = traj[k]
        vel = vels[k]
        # axis is the vector from microbe to target
        axis = target .- pos
        # project velocity onto axis
        v_drift[k] = dot(vel, axis) / norm(axis)
        if normalize
            v_drift[k] /= norm(vel)
        end # if
    end # for
    return v_drift
end # function


"""
    driftvelocity_direction(adf, target; normalize=false)
Evaluate the drift velocity of microbes along a direction `target`,
extracting their velocities from the agent dataframe `adf` (requires `:vel`
field).

Returns a matrix of instantaneous drift velocities with size
(nmicrobes, nsteps). By convention the drift velocity will be positive
for motion along the target direction.

If `normalize` is set to `true`, drift velocities are normalized by the
instantaneous speed of the microbe.
"""
function driftvelocity_direction(adf, target; normalize=false)
    vels = vectorize_adf_measurement(adf, :vel)
    if normalize
        v_drift = map(vel -> dot(vel, target) / norm(target) / norm(vel), vels)
    else
        v_drift = map(vel -> dot(vel, target) / norm(target), vels)
    end # if
    return v_drift
end # function