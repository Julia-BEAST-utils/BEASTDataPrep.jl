const CURRENT_NAME = [""]


function plot_transformed_data(df::DataFrame)

    nms = names(df)
    if nms[1] != TAXON_NAME
        error("first column of dataframe must be \"$(TAXON_NAME)\"")
    end
    n, p = size(df)

    plots = Matrix{Plot}(undef, p - 1, 3)



    for i = 2:p
        col = df[!, i]
        if !(eltype(col) <: Float64)
            error("all columns (except the first) must have elements of type Float64")
        end

        name = nms[i]
        CURRENT_NAME[1] = name

        p_original = plot(x=col,
                        Geom.histogram,
                        Guide.title(plot_title(i - 1, name, "untransformed")))
        plots[i - 1, 1] = p_original

        log_title = plot_title(i - 1, name, "log")
        logit_title = plot_title(i - 1, name, "logit")


        if can_transform_log(col)
            plots[i - 1, 2] = plot(x = log_transform(col),
                                Geom.histogram, Guide.title(log_title))
        else
            plots[i - 1, 2] = blank_plot(log_title)
        end
        if can_transform_logit(col)
            plots[i - 1, 3] = plot(x = logit_transform(col),
                                Geom.histogram, Guide.title(logit_title))
        else
            plots[i - 1, 3] = blank_plot(logit_title)
        end

    end

    plt = gridstack(plots)

    return plt
end

function plot_title(i::Int, name::String, transform::String)
    return "$i -> $name $transform"
end


function blank_plot(title::String)
    return plot(x=Float64[], Geom.histogram, Guide.title(title))
end

function can_transform_log(x::Vector{Float64})
    neg_ind = findfirst(y -> y < 0.0, x)
    return isnothing(neg_ind)
end


function log_transform(x::Vector{Float64})
    n = length(x)
    pos_inds = findall(y -> y > 0.0, x)
    log_x = zeros(n)

    if length(pos_inds) != n
        min = minimum(@view x[pos_inds])
        half_min = 0.5 * min
        log_half_min = log(half_min)
        zero_inds = setdiff(1:n, pos_inds)
        log_x[zero_inds] .= log_half_min

        @warn "0's in log-transform of $(CURRENT_NAME[1])." *
                " Converting all 0's to $half_min (half of the minium " *
                "value $min) prior to log-transforming."
    end

    for i in pos_inds
        log_x[i] = log(x[i])
    end
    return log_x
end

function can_transform_logit(x::Vector{Float64})
    bad_ind = findfirst(y -> y > 1.0 || y < 0.0, x)
    return isnothing(bad_ind)
end

function logit_transform(x::Vector{Float64})
    n = length(x)
    good_inds = findall(y -> y > 0.0 && y < 1.0, x)

    logit_x = zeros(n)

    if length(good_inds) != n
        good_x = @view x[good_inds]
        half_min = 0.5 * minimum(good_x)
        logit_half_min = logit(half_min)

        half_max = 0.5 * (1.0 + maximum(good_x))
        logit_half_max = logit(half_max)

        zero_inds = findall(isequal(0.0), x)
        logit_x[zero_inds] .= logit_half_min

        one_inds = findall(isequal(1.0), x)
        logit_x[one_inds] .= logit_half_max

        @warn "0's and/or 1's in logit-transform of $(CURRENT_NAME[1])." *
                " Converting all 0's to $half_min (half of the non-zero minium " *
                "value $min) and all 1's to $half_max (halfways between the " *
                "non-one maximum $half_max and 1.0) prior to logit-transforming."
    end

    for i in good_inds
        logit_x[i] = logit(x[i])
    end
    return logit_x
end

function logit(x::Float64)
    return log(x / (1 - x))
end