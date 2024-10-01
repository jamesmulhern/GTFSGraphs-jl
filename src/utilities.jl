function get_avg_transfer(arrivals::Vector{Int64}, departures::Vector{Int64})
    count = 0
    avg_headway = 0
    for arrival in arrivals
        idx = findfirst(departures .>= arrival)
        if !isnothing(idx)
            avg_headway += departures[idx] - arrival
            count += 1
        end
    end
    return avg_headway/count
end

function transfer_headway(src_df::SubDataFrame, dst_df::SubDataFrame)
    src_a = sort(duration.(view(src_df, :, :arrival_time)))
    src_d = sort(duration.(view(src_df, :, :departure_time)))
    dst_a = sort(duration.(view(dst_df, :, :arrival_time)))
    dst_d = sort(duration.(view(dst_df, :, :departure_time)))

    transfer = get_avg_transfer(src_a, dst_d)
    rev_transfer = get_avg_transfer(dst_a, src_d)

    return transfer, rev_transfer
end

function vertex_set_filters(gg::GTFSGraph{L,V,E}, label_vec::Vector{L}) where {L,V,E}
    
    #Create dict of sets for each field in label type
    prop_names = gg.vertex_label_params
    set_dict = Dict{Symbol, Set{String}}(prop_names .=> [Set{String}() for _ in prop_names])

    # Add labels to sets
    for label in label_vec
        for prop in prop_names
            push!(set_dict[prop],label[prop])
        end
    end

    # Create Filters
    filters = FilterDict()
    for p in prop_names
        merge!(filters, FilterDict(p => x -> in.(x, [set_dict[p]])))
    end

    return filters
end

function filter_vertexs(gg::GTFSGraph{L,V,E}, filter_func::Function) where {L,V,E}
    # Get list of all valid labels
    label_vec = Vector{L}()
    for label in labels(gg.g)
        if filter_func(gg.g[label])
            push!(label_vec, label)
        end
    end

    # Check that some labels are return
    if isempty(label_vec)
        warn("Vertex filter function $(filter_func) eliminated all vertexs.")
        return nothing
    end

    return label_vec
end


function make_tables(gg::GTFSGraph{L,V,E}) where {L,V,E}
    
    # Create vertex Table
    names = collect(fieldnames(V))
    pushfirst!(names,:label)
    nodes = DataFrame([name => [] for name in names])

    for label in labels(gg.g)
        d = Dict(n=>getfield(gg.g[label], n) for n in fieldnames(V))
        d = merge(d,Dict(:label=>string(label)))
        push!(nodes,d)
    end
    sort!(nodes, [:direction_id])

    # Create edge Table
    #dge = collect(edge_labels(g))[1]
    names = collect(fieldnames(E))
    pushfirst!(names,:dst)
    pushfirst!(names,:src)
    edges = DataFrame([name => [] for name in names])

    for (src, dst) in edge_labels(gg.g)
        d = Dict(n => getfield(gg.g[src,dst], n) for n in fieldnames(E)) 
        d = merge(d,Dict(:src=>string(src),:dst=>string(dst)))
        push!(edges,d)
    end
    sort!(edges, [:edge_type])

    return nodes, edges
end

function replace_named_tuples!(df::DataFrame)

    types = eltype.(eachcol(df))
    rmv_col = Vector{Symbol}()
    for (i, col) in enumerate(propertynames(df))
        #println("$col, $(types[i]), $(types[i] <: NamedTuple)")
        if col in [:stop_coord, :position]
            fields = fieldnames(typeof(df[1,col]))
            insertcols!(df, col, [f => [r[col][f] for r in eachrow(df)] for f in fields]...)
            push!(rmv_col, col)
        end
    end

    # Remove orginal columns
    select!(df, Not(rmv_col)) 
end


function save_graph(base_name::String, gg::GTFSGraph, vertex_params::Vector{Symbol}=Vector{Symbol}(), edge_params::Vector{Symbol}=Vector{Symbol}())

    nodes, edges = make_tables(gg)
    n = nrow(nodes)
    e = nrow(edges)
    
    insertcols!(nodes, propertynames(nodes)[1], :index => collect(1:n))
    insertcols!(edges, propertynames(edges)[1], :index => collect(1:e))

    replace!.([edges.src], nodes.label .=> string.(nodes.index))
    replace!.([edges.dst], nodes.label .=> string.(nodes.index))

    replace_named_tuples!(nodes)

    # Add filters for columns

    n_file = "$(base_name).nodes.csv"
    e_file = "$(base_name).edges.csv"

    CSV.write(n_file, nodes)
    CSV.write(e_file, edges)
end