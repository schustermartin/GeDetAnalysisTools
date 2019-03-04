
mutable struct Detector
  name::AbstractString
  n_channels::Int
  channel_display_order::Array{Int}  # channel display order
  channel_display_layout::Array{Int} # channel display layout
  chn_idx_neighbour_left::Array
  chn_idx_neighbour_right::Array
  chn_idx_neighbour_top::Array
  chn_idx_neighbour_bottom::Array
  crystal_radius::Float64
  crystal_length::Float64
  channel_plot_layout::Function

  function Detector()
    new("noname", 2, [1,2], [1,2], [], [], [], [], 0, 0, identity)
  end
end

function SIEGFRIED_3()
    s3 = Detector()
    s3.name = "Siegfried 3"
    s3.n_channels = 19
    s3.crystal_length = 70
    s3.crystal_radius = 37.5

    s3.channel_display_order = [1,7,8,9,13,14,15,19,20,21,22,23,24,16,17,18,10,11,12]
    # -> for @layout [  chn1{0.5w} chn20{0.5w}
            # chn2 chn3 chn4 chn5 chn6 chn7
            # chn8 chn9 chn10 chn11 chn12 chn13
            # chn14 chn15 chn16 chn17 chn18 chn19]

    s3.channel_display_layout = [4,6]
    s3.chn_idx_neighbour_left =   [ 18, 1, 2, 15, 4, 5, 12, 7, 8, 3, 16, 17, 6, 13, 14, 9, 10, 11 ] .+ 1
    s3.chn_idx_neighbour_right =  [ 2, 3, 16, 5, 6, 13, 8, 9, 10, 11, 12, 7, 14, 15, 4, 17, 18, 1] .+ 1
    s3.chn_idx_neighbour_top =    [ -1, -1, -1, 1, 2, 3, 4, 5, 6, 13, 14, 15, 16, 17, 18, -1, -1, -1 ] .+ 1
    s3.chn_idx_neighbour_bottom = [ 4, 5, 6, 7, 8, 9, -1, -1, -1, -1, -1, -1, 10, 11, 12, 13, 14, 15 ] .+ 1
    return s3
end

function SUSIE()
    susie = Detector()
    susie.name = "Susie"
    susie.n_channels = 20
    susie.crystal_length = 70
    susie.crystal_radius = 37.5
    # susie.channel_display_order = [1,2,6,10,3,7,11,4,8,12,16,20,24,15,19,23,14,18,22] old ROOT
    susie.channel_display_order = [1, 20, 19, 18, 14, 13, 12, 8, 7, 6, 5, 4, 3, 11, 10, 9, 17, 16, 15, 2] .+ 1
    susie.channel_display_order[1] = 1
    susie.channel_display_order[20] = 2
    susie.channel_display_layout = [4, 6]
    susie.chn_idx_neighbour_left =   [ missing, 3, 4, 17, 6, 7, 14, 9, 10, 11, 12, 13, 8, 15, 16, 5, 18, 19, 1, missing ] 
    susie.chn_idx_neighbour_right =  [ missing, 19, 2, 3, 16, 5, 6, 13, 8, 9, 10, 11, 12, 7, 14, 15, 4, 17, 18, missing ] 
    susie.chn_idx_neighbour_top =    [ missing, 5, 6, 7, 8, 9, 10, 20, 20, 20, 20, 20, 20, 11, 12, 13, 14, 15, 16, missing] 
    susie.chn_idx_neighbour_bottom = [ missing, missing, missing, missing, 2, 3, 4, 5, 6, 7, 14, 15, 16, 17, 18, 19, missing, missing, missing, [8, 9, 10, 11, 12, 13]] 

    function layout_func()
      return @layout [ chn1 chn20 info{0.66666666667w}
                       grid(3, 6){0.75h} ]
                        # chn2 chn3 chn4 chn5 chn6 chn7
                        # chn8 chn9 chn10 chn11 chn12 chn13
                        # chn14 chn15 chn16 chn17 chn18 chn19]
    end
    susie.channel_plot_layout = layout_func

    return susie
end


function segBEGe()
    s3 = Detector()
    s3.name = "segBEGe"
    s3.n_channels = 5
    s3.channel_display_order = [1,4,5,6,3]
    # s3.channel_display_order = [1,2,3,4,5]
    s3.channel_display_layout = [3,2]
    function layout_func()
        return @layout grid(3,2)
    end
    s3.channel_plot_layout = layout_func
    return s3
end
