function circle(x, y, r=1; n=30)
    θ = 0:360÷n:360
    Plots.Shape(r*sind.(θ) .+ x, r*cosd.(θ) .+ y)
end

@recipe function f(m::AbstractMicrobe{2})
    x₀ = m.pos[1]
    y₀ = m.pos[2]
    r = m.radius == 0 ? 5 : m.radius * 10
    dx, dy = m.vel ./ (norm(m.vel)) .* r

    label --> false
    
    # draw circle for body
    @series begin
        fillcolor --> 1
        circle(x₀, y₀, r)
    end
    # draw arrow for velocity
    @series begin
        label := false
        seriestype := :line
        arrows := true
        linecolor --> 2
        linewidth --> 2 
        [x₀, x₀+dx], [y₀, y₀+dy]
    end
end # recipe

@recipe function f(microbes::T) where {S<:AbstractMicrobe{2},T<:AbstractVector{S}}
    for m in microbes
        x₀ = m.pos[1]
        y₀ = m.pos[2]
        r = m.radius == 0 ? 5 : m.radius * 10
        dx, dy = m.vel ./ (norm(m.vel)) .* r

        label --> false
        
        # draw circle for body
        @series begin
            fillcolor --> 1
            circle(x₀, y₀, r)
        end
        # draw arrow for velocity
        @series begin
            label := false
            seriestype := :line
            arrows := true
            linecolor --> 2
            linewidth --> 2
            [x₀, x₀+dx], [y₀, y₀+dy]
        end
    end # for
end # recipe