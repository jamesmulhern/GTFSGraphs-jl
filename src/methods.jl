function connect_vertexs!(gg::GTFSGraph{L,V,E}, radius::Float64, filter_func::Function=Returns(true)) where {L, V, E}
    # Get list of all valid labels
    valid_labels = filter_vertexs(gg, filter_func)

    # Compute nearby vertexs
    positions = hcat(collect.([gg.g[label].stop_coord for label in valid_labels])...)
    tree = BallTree(positions, gg.distance_function)
    near_idx_list = inrange(tree, positions, radius)

    for (i, near_idxs) in pairs(near_idx_list)
        # Loop over local vector if near indexs
        for idx in near_idxs
            #Skip if referencing itself
            i == idx ? continue : nothing

            src = valid_labels[i]
            dst = valid_labels[idx]

            # Check if key already exists
            haskey(gg.g, src, dst) ? continue : nothing

            # Compute headway 
            distance = gg.distance_function(gg.g[src].stop_coord, gg.g[dst].stop_coord)

            # Add forward edge and reverse edge
            gg.g[src,dst] = E()
            gg.g[src,dst].shape_distance = distance
            gg.g[src,dst].edge_type = "transfer"

            gg.g[dst,src] = E()
            gg.g[dst,src].shape_distance = distance
            gg.g[dst,src].edge_type = "transfer"
        end
    end
    return nothing
end

function _set_transfer_weight!(e::AbstractEdgeData; kwargs...)
    src_df = kwargs[:src_df]
    dst_df = kwargs[:dst_df]

    fwd_transfer, _ = transfer_headway(src_df, dst_df)
    setfield!(e, :s2s_time, fwd_transfer)
    return nothing
end

function link_vertexs2!(gg::GTFSGraph{L,V,E}, radius::Float64, filter_func::Function=Returns(true)) where {L, V, E}
    connect_vertexs!(gg, radius, filter_func)
    new_func = x -> filter_func(x) && ==(x.edge_type, "transfer")
    apply_edge_function!(gg, [:arrival_time, :departure_time], _set_transfer_weight!, new_func)

end

function link_vertexs!(gg::GTFSGraph{L,V,E}, radius::Float64, filter_func::Function=Returns(true)) where {L, V, E}

    # Get list of all valid labels
    valid_labels = Vector{L}()
    for label in labels(gg.g)
        if filter_func(gg.g[label])
            push!(valid_labels, label)
        end
    end

    # Check that some labels are return
    if isempty(valid_labels)
        warn("Vertex filter function $(filter_func) eliminated all vertexs. No links were added")
        return nothing
    end

    params = [:arrival_time, :departure_time, gg.vertex_label_params...]
    filters = vertex_set_filters(gg, valid_labels)
    df = query(gg.q, params, filters)
    gdf = groupby(df, gg.vertex_label_params)
    #gdf = query_stop_times(gg, valid_labels)

    # Compute nearby vertexs
    positions = hcat(collect.([gg.g[label].stop_coord for label in valid_labels])...)
    tree = BallTree(positions, gg.distance_function)
    near_idx_list = inrange(tree, positions, radius)
    
    # Loop over all points
    for (i, near_idxs) in pairs(near_idx_list)

        # Loop over local vector if near indexs
        for idx in near_idxs
            #Skip if referencing itself
            if i == idx
                continue
            end

            src = valid_labels[i]
            dst = valid_labels[idx]

            # Check if key already exists
            if haskey(gg.g, src, dst)
                continue
            end

            # Check df has key  # TODO maybe remove?
            if !haskey(gdf, dst)
                continue
            end

            # Get Data from grouped DataFrame
            src_df = view(gdf[src], :, [:arrival_time,:departure_time])
            dst_df = view(gdf[dst], :, [:arrival_time,:departure_time])

            # Compute headway 
            fwd_transfer, rev_transfer = transfer_headway(src_df, dst_df)
            distance = gg.distance_function(gg.g[src].stop_coord, gg.g[dst].stop_coord)

            # Add forward edge and reverse edge
            gg.g[src,dst] = E()
            gg.g[src,dst].shape_distance = distance
            gg.g[src,dst].s2s_time = fwd_transfer
            gg.g[src,dst].edge_type = "transfer"

            gg.g[dst,src] = E()
            gg.g[dst,src].shape_distance = distance
            gg.g[dst,src].s2s_time = rev_transfer
            gg.g[dst,src].edge_type = "transfer"
        end
    end
    return nothing
end


function _add_headway_and_dwell!(vertex_data::AbstractVertexData; kwargs...)
    df = kwargs[:df]
    arrive = sort(duration.(view(df, :, :arrival_time)))
    depart = sort(duration.(view(df, :, :departure_time)))

    avg_dwell = get_avg_transfer(arrive, depart)
    #avg_headway = get_avg_transfer(depart[1:end-1], depart[2:end])
    avg_headway = (60*60*(20-6))/length(depart)

    setfield!(vertex_data, :dwell_time, avg_dwell)
    setfield!(vertex_data, :headway, avg_headway)
    return nothing
end

function add_headway_and_dwell!(gg::GTFSGraph{L, V, E}, vertex_filter::Function=Returns(true)) where {L,V,E}

    apply_vertex_function!(gg, [:arrival_time, :departure_time], _add_headway_and_dwell!, vertex_filter)
    return nothing
end

function _add_running_time!(edge_data::AbstractEdgeData; kwargs...)

    src_df = kwargs[:src_df]
    dst_df = kwargs[:dst_df]

    src_d = sort(duration.(view(src_df, :, :departure_time)))
    dst_a = sort(duration.(view(dst_df, :, :arrival_time)))

    avg_time = get_avg_transfer(src_d, dst_a)

    # Update Edge weight
    setfield!(edge_data, :running_time, avg_time)
    return nothing
end

function add_running_time!(gg::GTFSGraph{L, V, E}, edge_filter::Function=Returns(true)) where {L, V, E}
    apply_edge_function!(gg, [:arrival_time, :departure_time], _add_running_time!, edge_filter)
    return nothing
end

function _add_s2s_time!(edge_data::AbstractEdgeData; kwargs...)
    setfield!(edge_data, :s2s_time, kwargs[:src_data].dwell_time + edge_data.running_time)
    return nothing
end

function add_s2s_time!(gg::GTFSGraph{L,V,E}, edge_filter::Function=Returns(true)) where {L,V,E}
    apply_edge_function!(gg, Vector{Symbol}(), _add_s2s_time!, edge_filter)
    return nothing
end


    # labels = Vector{L}
    # for (src, dst) in edge_labels(gg.g)
    #     if edge_filter(gg.g[src,dst])
    #         push!(labels, src)
    #         push!(labels, dst)
    #     end
    # end

    # # Query data related to stops
    # params = [:arrival_time, :departure_time, gg.vertex_label_params...]
    # filters = vertex_set_filters(labels)
    # df = query(gg.q, params, filters) #Need to apply date filters
    # gdf = groupby(df, gg.vertex_label_params)


    # # Loop over edges
    # for (src, dst) in edge_labels(gg.g) 
    #     #Skip edges not in dataset
    #     if !haskey(gdf, src) || !haskey(gdf, dst)
    #         continue
    #     end

    #     src_d = sort(duration.(view(gdf[src], :, :departure_time)))
    #     dst_a = sort(duration.(view(gdf[dst], :, :arrival_time)))

    #     avg_time = get_avg_transfer(src_d, dst_a)

    #     # Update Edge weight
    #     setfield!(gg.g[src,dst], :weight, avg_time)
    # end