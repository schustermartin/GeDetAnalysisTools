dictvalues(d::Dict) = (collect(values(d))...,)
dictkeys(d::Dict) = (collect(keys(d))...,)

namedtuple(d::Dict{String,T}) where {T} =
    NamedTuple{Symbol.(dictkeys(d))}(dictvalues(d))

function read_parameters()
    Filename = dirname(pathof(GeDetAnalysisTools))*"/filters/parameters.json"
    Dict     = JSON.parse(String(read(Filename)))
    Keys     = String.(keys(Dict))
    i = 1
    while i <= length(Keys)
        Dict[Keys[i]] = namedtuple(Dict[Keys[i]])
        i += 1
    end
    return sort(Dict)
end

function add_parameter(Item_Name::String, p::NamedTuple)
    Filename = dirname(pathof(GeDetAnalysisTools))*"/filters/parameters.json"
    D        = read_parameters(Filename)
    D[Item_Name] = p
    JSON_String  = json(sort(D))
    open(Filename, "w") do f
        JSON.print(f, D, 4) end
    return sort(D)
end

edit_parameter(Item_Name::String, p::NamedTuple) = add_parameter(Item_Name::String, p::NamedTuple)

function rm_parameter(Item_Name::String)
    Filename = dirname(pathof(GeDetAnalysisTools))*"/filters/parameters.json"
    D = read_parameters(Filename)
    delete!(D, Item_Name)
    JSON_String = json(sort(D))
    open(Filename, "w") do f
        JSON.print(f, D, 4) end
    return sort(D)
end