
module ComposeAtoms

using Compose, Colors, JuLIP

import Base.display
import Compose: compose, context

export compose, compose_layers, display_layers


plotschemes = Dict(
   :default =>
   Dict( :atcol    => "tomato",
         :qmcol    => "cyan",
         :bufcol   => "purple",
         :bondcol  => "darkblue",
         :bdrycol => "darkgrey",
         :atradius => 0.2,
         :bdryradius => 0.1,
         :atborder => 0.3,
         :bondwidth => 0.3,
         :axisbuffer => 0.05,
         :atcol2 => "aquamarine2"
         ),

   :presentation =>
      Dict( :atcol    => RGB(78/255,115/255,174/255),    # blue
            :qmcol    => RGB(194/255, 79/255, 84/255),   # red
            :bufcol   => "purple",   # RGB(129/255,116/255,176/255),   # purple
            :bondcol  => RGB(0.35, 0.35, 0.35),     # grey
            :bdrycol => RGB(0.7, 0.7, 0.7),     # grey
            :atradius => 0.2,
            :bdryradius => 0.1,
            :atborder => 0.3,
            :bondwidth => 0.3,
            :axisbuffer => 0.05
            ),

   :P2 =>
      Dict( :atcol    => RGB(100/255,130/255,220/255),    # blue
            :qmcol    => RGB(220/255, 120/255, 120/255),   # red
            :bufcol   => RGB(170/255,150/255,210/255),   # RGB(129/255,116/255,176/255),   # purple
            :bondcol  => RGB(0.35, 0.35, 0.35),     # grey
            :bdrycol => RGB(0.7, 0.7, 0.7),     # grey
            :atradius => 0.2,
            :bdryradius => 0.1,
            :atborder => 0.3,
            :bondwidth => 0.3,
            :axisbuffer => 0.05
            ),

   :layers =>
      Dict( :atcol    => "tomato",
            :qmcol    => "cyan",
            :bufcol   => "purple",
            :bondcol  => RGB(0.35, 0.35, 0.35),
            :bdrycol => "darkgrey",
            :atradius => 0.2,
            :atborder => 0.3,
            :bondwidth => 0.3,
            :axisbuffer => 0.05,
            :atcolors => ["tomato", "aquamarine3", "purple"],
            :atradii => [0.4, 0.35, 0.3],
            :bdryradius => 0.1
   ),
)


plotscheme(s) = plotschemes[s]

function analyse_kwargs!(scheme; kwargs...)
   for (key, val) in kwargs
      scheme[key] = val
   end
end

coordinates(X::Matrix) = X[1,:][:], X[2,:][:]
coordinates(X::Vector{JVecF}) = coordinates(mat(X))

"return a reasonable choice for the plotting axis"
function autoaxis(X, scheme)
   x, y = coordinates(X)
    xLim = [extrema(x)...]
    width = xLim[2]-xLim[1]
    yLim = [extrema(y)...]
    height = yLim[2]-yLim[1]
    buffer_factor = scheme[:axisbuffer]
    buffer = min( 1.0, 0.05 * width, 0.05 * height )
    xLim[1] -= buffer; xLim[2] += buffer
    yLim[1] -= buffer; yLim[2] += buffer
    return [xLim[1], xLim[2], yLim[1], yLim[2]]
end

"turn axis into UnitBox for creating contexts"
_ub_(axis) = UnitBox(axis[1], axis[3], axis[2]-axis[1], axis[4]-axis[3])

context(axis::Vector) = context(units=_ub_(axis))

relative_height(ax, width) = (ax[4]-ax[3])/(ax[2]-ax[1]) * width
relative_width(ax, height) = (ax[2]-ax[1])/(ax[4]-ax[3]) * height


function compose_atoms(X, axis, r, atcol, atborder, bordercol)
   # return an empty context if we have an empty set of atoms
   if length(X) == 0
      return (context(axis),)
   end
   x, y = coordinates(X)
   r = ones(length(x)) * r
   return ( context(axis), circle(x, y, r),
              stroke(bordercol), linewidth(atborder), fill(atcol) )
end


compose_bg(axis, bg) = (
   bg ? (compose(context(axis), rectangle()), fill("white")) :
   (context(axis),) )


function unique_bonds(at::AbstractAtoms, rcut; X = nothing)
   bc = pbc(at)
   set_pbc!(at, false)
   if X != nothing
      X0 = positions(at)
      set_positions!(at, X)
   end
   nlist = neighbourlist(at, rcut)
   ij = sort([nlist.i nlist.j], 2)
   ij = unique(ij, 1)
   if X != nothing
      set_positions!(at, X0)
   end
   set_pbc!(at, bc)
   return (ij[:,1], ij[:, 2])
end



# edges2path(X, B) = insert_nans( [X[:, B[1]]; X[:, B[2]]] )
# insert_nans(P) = reshape( [P; NaN * ones(2, size(P, 2))], 2, 3 * size(P, 2) )

function compose_bonds(X, B, scheme, axis)
   i, j = B[1], B[2]
   x, y = coordinates(X)
   P = NTuple{2, Float64}[]
   for n = 1:length(i)
      append!(P, [ (x[i[n]], y[i[n]]), (x[j[n]], y[j[n]]), (NaN, NaN)])
   end
   # P = mat2points(reshape(P, 2, length(P) ÷ 2))
   # return the required tuple to compose
   return ( context(axis), line(P), stroke(scheme[:bondcol]), linewidth(scheme[:bondwidth]) )
end



function compose(at::AbstractAtoms;
                  X=nothing, scheme = plotscheme[:default],
                  axis = autoaxis(X, scheme),
                  rcut = rnn(at) * 1.2, bg = false,
                  Iqm = Int[], Ibdry = Int[], Ibuf = Int[],
                  kwargs... )

   analyse_kwargs!(scheme; kwargs...)
   i, j = unique_bonds(at, rcut; X=X)
   if X == nothing
      X = positions(at)
   end

   Idom = collect(1:length(at))
   Imm = setdiff(Idom, union(Iqm, Ibdry, Ibuf))

   ctx = compose( context(axis),
      compose_atoms( X[Imm], axis, scheme[:atradius], scheme[:atcol], scheme[:atborder], scheme[:bondcol] ),
      compose_atoms( X[Iqm], axis, scheme[:atradius], scheme[:qmcol], scheme[:atborder], scheme[:bondcol] ),
      compose_atoms( X[Ibuf], axis, scheme[:atradius], scheme[:atcol], scheme[:atborder], scheme[:qmcol] ),
      compose_atoms( X[Ibdry], axis, scheme[:bdryradius], scheme[:bdrycol], scheme[:atborder], scheme[:bondcol] ),
      compose_bonds( X, (i, j), scheme, axis )
                  # compose_bg(axis, bg)
               )
   return ctx
end



function compose_layers( at::AbstractAtoms, zlayers;
                        X=nothing, scheme = plotscheme[:layers],
                        axis = autoaxis(X, scheme),
                        rcut = rnn(at) * 1.2, bg = false,
                        Iqm = Int[], Ibdry = Int[], Ibuf = Int[],   # ignore all except Ibrdy for now
                        kwargs... )

   analyse_kwargs!(scheme; kwargs...)
   i, j = unique_bonds(at, rcut; X=X)
   if X == nothing
      X = positions(at)
   end

   # add atoms to one of the two layers
   z = mat(X)[3,:]
   II = [ Int[] for j = 1:length(zlayers) ]
   for j = 1:length(z)
      _, i = findmin(abs(z[j] - zlayers))
      push!(II[i], j)
   end
   for Ij in II
      Ij = setdiff(Ij, Ibdry)
   end

   cat = [ compose_atoms( X[II[j]], axis, scheme[:atradii][j], scheme[:atcolors][j],
                          scheme[:atborder], scheme[:bondcol] )  for j = 1:length(II) ]

   ctx = compose( context(axis), cat...,
      compose_atoms( X[Ibdry], axis, scheme[:bdryradius], scheme[:bdrycol], scheme[:atborder], scheme[:bondcol] ),
      compose_bonds( X, (i, j), scheme, axis )
               )
   return ctx
end



function _display(ctx, filename, display, plotwidth, plotheight, Img)
   if (filename == nothing) || (display == true)
      draw( SVG(plotwidth, plotheight), ctx )
   end
   if filename != nothing
      draw( Img(filename, plotwidth, plotheight), ctx )
   end
end


"""
`display(at::AbstractAtoms; kwargs...)`

## Keyword arguments
* X : positions (default, reference positions)
* axis : default is an autoaxis
* plotwidth: default is 12cm
* Img : image type, default is SVG
* rnn : bond-length, default given by `rnn(at)` (cf JuLIP)
"""
display(at::AbstractAtoms;
               X=positions(at), scheme = plotscheme(:default),
               axis = autoaxis(X, scheme), display = true,
               plotwidth = 12cm, plotheight= relative_height(axis, plotwidth),
               Img=SVG, filename = nothing, kwargs... ) =
   _display( compose(at; X=X, axis=axis, scheme=scheme, kwargs...),
             filename, display, plotwidth, plotheight, Img )


display_layers(at::AbstractAtoms, zlayers;
                  X=positions(at), scheme = plotscheme(:layers),
                  axis = autoaxis(X, scheme), display = true,
                  plotwidth = 12cm, plotheight= relative_height(axis, plotwidth),
                  Img=SVG, filename = nothing, kwargs... ) =
   _display( compose_layers(at, zlayers; X=X, axis=axis, scheme=scheme, kwargs...),
            filename, display, plotwidth, plotheight, Img )



end
