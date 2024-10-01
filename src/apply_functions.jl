# Funcitons that allow the user to apply a custom funciton to edges or vertexs

#
# Vertex Functions
#

function apply_vertex_function!(gg::GTFSGraph{L,V,E}, params::Vector{Symbol}, func::Function, label_vec::Vector{L}) where {L, V, E}
    
    # TODO Check vertex function is of correct type
    # my_ertex_func!(vertex_data::V; data::SubDataFrame, label::L)
    # inputs - VertexData::V - required
    #          query_results::SubDataFrame
    #          label=L -optional
    # return - nothing, mutate input data
    
    # Query data    
    all_params = vcat(params, gg.vertex_label_params) |> unique
    filters = vertex_set_filters(gg, label_vec)
    df = query(gg.q, all_params, filters)
    gdf = groupby(df, gg.vertex_label_params)
    

    # Apply function
    for label in label_vec
        func(gg.g[label]; label=label, df=gdf[label])
    end
    return nothing
end

function apply_vertex_function!(gg::GTFSGraph{L,V,E}, params::Vector{Symbol}, func::Function, filter_func::Function=Returns(true)) where {L,V, E}

    # Create list of valid labels
    label_vec = filter_vertexs(gg, filter_func)
    # for label in labels(gg.g)
    #     if filter_func(gg.g[label])
    #         push!(label_vec, label)
    #     end
    # end

    # # warn and return if no labels are
    # if isempty(label_vec) 
    #     warn("Filtering function $filter_func eliminated all vertexs, vertex function was not applied to any vertexs")
    #     return nothing
    # end

    apply_vertex_function!(gg, params, func, label_vec)
    return nothing
end

#
# Edge Functions
#

function apply_edge_function!(gg::GTFSGraph{L,V,E}, params::Vector{Symbol}, func::Function, edge_pairs::Vector{Pair{L,L}}) where {L,V,E}

    # TODO check function layout
    # my_edge_func(EdgeData::E; src_df::SubDataFrame, dst_df:SubDataFrame, src::L, dst::L, )

    
    # Get vector of labels
    label_vec = vcat(getindex.(edge_pairs,1), getindex.(edge_pairs,2)) |> unique

    # Perfrom query
    all_params = vcat(params, gg.vertex_label_params) |> unique
    filters = vertex_set_filters(gg, label_vec)
    df = query(gg.q, all_params, filters)
    gdf = groupby(df, gg.vertex_label_params)

    for (src,dst) in edge_pairs
        func(gg.g[src,dst]; src=src, src_df=gdf[src], src_data=gg.g[src], dst=dst, dst_df=gdf[dst], dst_data=gg.g[dst])
    end
    return nothing
end

function apply_edge_function!(gg::GTFSGraph{L,V,E}, params::Vector{Symbol}, func::Function, filter_func::Function=Returns(true)) where {L,V,E}
    edge_pairs = Vector{Pair{L,L}}()

    for (src,dst) in edge_labels(gg.g)
        if filter_func(gg.g[src,dst])
            push!(edge_pairs, src=>dst)
        end
    end

    # warn and return if no labels are
    if isempty(edge_pairs) 
        warn("Filtering function $filter_func eliminated all edges, edge function was not applied to any edges")
        return nothing
    end

    apply_edge_function!(gg, params, func, edge_pairs)
    return nothing
end