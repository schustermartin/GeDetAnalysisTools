export Detector

mutable struct Detector
  name::AbstractString
  n_channels::Int
  channel_display_order::Array{Int}  # channel display order
  channel_display_layout::Array{Int} # channel display layout
  chn_idx_neighbour_left::Array{Int}
  chn_idx_neighbour_right::Array{Int}
  chn_idx_neighbour_top::Array{Int}
  chn_idx_neighbour_bottom::Array{Int}

  function Detector()
    new("noname", 2, [1,2], [1,2], [], [], [], [])
  end
end

export SIEGFRIED_3
function SIEGFRIED_3()
    s3 = Detector()
    s3.name = "Siegfried 3"
    s3.n_channels = 19

    s3.channel_display_order = [1,7,8,9,13,14,15,19,20,21,22,23,24,16,17,18,10,11,12]
    # -> for @layout [  chn1{0.5w} chn20{0.5w}
            # chn2 chn3 chn4 chn5 chn6 chn7
            # chn8 chn9 chn10 chn11 chn12 chn13
            # chn14 chn15 chn16 chn17 chn18 chn19]

    s3.channel_display_layout = [4,6]
    s3.chn_idx_neighbour_left = [ 18, 1, 2, 15, 4, 5, 12, 7, 8, 3, 16, 17, 6, 13, 14, 9, 10, 11 ] .+ 1
    s3.chn_idx_neighbour_right = [ 2, 3, 16, 5, 6, 13, 8, 9, 10, 11, 12, 7, 14, 15, 4, 17, 18, 1] .+ 1
    s3.chn_idx_neighbour_top = [ -1, -1, -1, 1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, -1, -1, -1 ] .+ 1
    s3.chn_idx_neighbour_bottom = [ 4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, 10, 11, 12, 13, 14, 15 ] .+ 1
    return s3
end

export SUSIE
function SUSIE()
    s3 = Detector()
    s3.name = "Susie"
    s3.n_channels = 20
    # s3.channel_display_order = [1,2,6,10,3,7,11,4,8,12,16,20,24,15,19,23,14,18,22] old ROOT
    s3.channel_display_order = [1, 3, 4, 5, 9, 10, 11, 15, 16, 17, 18, 19, 20, 12, 13, 14, 6, 7, 8,2]
    s3.channel_display_layout = [4,6]
    s3.chn_idx_neighbour_left = [ 18, 1, 2, 15, 4, 5, 12, 7, 8, 3, 16, 17, 6, 13, 14, 9, 10, 11, -1 ] .+ 1
    s3.chn_idx_neighbour_right = [ 2, 3, 16, 5, 6, 13, 8, 9, 10, 11, 12, 7, 14, 15, 4, 17, 18, 1, -1] .+ 1
    s3.chn_idx_neighbour_top = [ 19, 19, 19, 1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, 19, 19, 19, -1] .+ 1
    s3.chn_idx_neighbour_bottom = [ 4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, 10, 11, 12, 13, 14, 15, -1 ] .+ 1
    return s3
end


export segBEGe
function segBEGe()
    s3 = Detector()
    s3.name = "segBEGe"
    s3.n_channels = 5
    s3.channel_display_order = [1,4,5,6,3]
    s3.channel_display_layout = [3,2]
    return s3
end
