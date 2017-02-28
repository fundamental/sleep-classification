#assumes model is the random forest
function get_feature_use(model)
    feats = zeros(1000)
    for tree in model.trees
        get_feature_use_recur(feats, tree)
    end
    feats
end

function get_feature_use_recur(feats, node)
    if(depth(node) < 1)
        return
    end
    if(typeof(node) != DecisionTree.Leaf)
        feats[node.featid] += 1
        get_feature_use_recur(feats, node.left)
        get_feature_use_recur(feats, node.right)
    end
end

#assumes model is the random forest
function get_class_difficulty(model)
    feats = Any[]
    push!(feats, [])
    push!(feats, [])
    push!(feats, [])
    push!(feats, [])
    push!(feats, [])
    for tree in model.trees
        get_class_difficulty_recur(feats, tree, 0)
    end
    feats
end

function get_class_difficulty_recur(feats, node, depth)
    if(typeof(node) != DecisionTree.Leaf)
        get_class_difficulty_recur(feats, node.left,  depth+1)
        get_class_difficulty_recur(feats, node.right, depth+1)
    else
        push!(feats[node.majority], depth)
    end
end
