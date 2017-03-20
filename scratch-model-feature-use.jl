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

function apply_model_prob(feats, model)
    prob = zeros(5, size(feats,2))
    trees = length(model.trees)
    for i = 1:size(feats,2)
        for j=1:trees
            out  = apply_tree(feats[:,i], model.tree)
            prob[out, i] += 1.0/trees
        end
    end
    prob
end

