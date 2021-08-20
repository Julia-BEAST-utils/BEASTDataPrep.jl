const CURRENT_NAME = [""]

################################################################################
## Makes a plot with all relevant transforms to visualize the data
################################################################################

function plot_transformed_data(df::DataFrame;
                               svg_path::String = "checkTransforms.svg",
                               overwrite::Bool = false)

    if !overwrite && isfile(svg_path)
        throw(ArgumentError("file '$svg_path' already exists. Either specify " *
                    "'overwrite=true' or 'svg_path=<some other file name>'."))
    end

    nms = names(df)

    n, p = size(df)
    relevant_columns = findall(x -> eltype(x) <: Real, [df[!, i] for i = 1:p])

    p_rel = length(relevant_columns)
    plots = Matrix{Plot}(undef, p_rel, 3)

    p_df = DataFrame(trait = names(df)[relevant_columns],
            untransformed = fill(NaN, p_rel),
            log = fill(NaN, p_rel),
            logit = fill(NaN, p_rel))



    for row in 1:length(relevant_columns)
        i = relevant_columns[row]
        col = df[!, i]
        if !(eltype(col) <: Real)
            error("all columns (except the first) must have elements of type Real")
        end

        name = nms[i]
        CURRENT_NAME[1] = name

        col = col[findall(!isnan, col)]

        p = pvalue(JarqueBeraTest(col))
        p_df.untransformed[row] = p

        p_original = plot(x=col,
                        Geom.histogram,
                        Guide.title(plot_title(row, name, "untransformed")))# p = $p")))
        plots[row, 1] = p_original

        log_title = plot_title(row, name, "log")
        logit_title = plot_title(row, name, "logit")


        if can_transform_log(col)
            logged_col = log_transform(col)
            p = pvalue(JarqueBeraTest(logged_col))
            p_df.log[row] = p
            # log_title = log_title * " p = $p"
            plots[row, 2] = plot(x = logged_col,
                                Geom.histogram, Guide.title(log_title))
        else
            plots[row, 2] = blank_plot(log_title)
        end
        if can_transform_logit(col)
            logit_col = logit_transform(col)
            p = pvalue(JarqueBeraTest(logit_col))
            p_df.logit[row] = p
            # logit_title = logit_title * " p = $p"
            plots[row, 3] = plot(x = logit_col,
                                Geom.histogram, Guide.title(logit_title))
        else
            plots[row, 3] = blank_plot(logit_title)
        end

    end

    plt = gridstack(plots)
    draw(SVG(svg_path, 10inch, p * 2inch), plt)
    return p_df
end

function plot_title(i::Int, name::String, transform::String)
    return "$i -> $name $transform"
end


function blank_plot(title::String)
    return plot(x=Float64[], Geom.histogram, Guide.title(title))
end

################################################################################
## General tools for transforming data
################################################################################

function can_transform_log(x::Vector{<:Real})
    neg_ind = findfirst(y -> y < 0.0, x)
    return isnothing(neg_ind)
end


function log_transform(x::Vector{<:Real})
    n = length(x)
    pos_inds = findall(y -> y > 0.0, x)
    log_x = fill(NaN, n)

    if length(pos_inds) != n
        min = minimum(@view x[pos_inds])
        half_min = 0.5 * min
        log_half_min = log(half_min)
        zero_inds = findall(isequal(0.0), x)
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

function can_transform_logit(x::Vector{<:Real})
    bad_ind = findfirst(y -> y > 1.0 || y < 0.0, x)
    return isnothing(bad_ind)
end

function logit_transform(x::Vector{<:Real})
    n = length(x)
    good_inds = findall(y -> y > 0.0 && y < 1.0, x)

    logit_x = fill(NaN, n)

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

################################################################################
## transforming the dataframe
################################################################################

function transform_data!(df::DataFrame, transform::Function;
                        inds::AbstractVector{Int} = 1:(size(df, 2) - 1))
    nms = names(df)
    for ind in inds
        CURRENT_NAME[1] = nms[ind + 1]
        col = df[!, ind + 1]
        col .= transform(col)
    end
end


function log_transform_data!(df::DataFrame;
                             inds::AbstractVector{Int} = 1:(size(df, 2) - 1))
    transform_data!(df, log_transform, inds = inds)
end

function standardize_data!(df::DataFrame)
    for i = 2:size(df, 2)
        col = df[findall(!isnan, df[!, i]), i]
        μ = mean(col)
        σ = std(col, mean=μ, corrected=false)
        df[!, i] .-= μ
        df[!, i] ./= σ
    end
end

