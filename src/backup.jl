# function compose_core(X, B, rdef)
#    nX = size(X,2)
#    r = sumabs2(X, 1) |> sqrt
#    nneigs = zeros(Int, nX)
#    for n = 1:length(B[1])
#       nneigs[B[1][n]] += 1
#       nneigs[B[2][n]] += 1
#    end
#    Icore = find( (nneigs .!= 6) .* (r .<= rdef) )
# end


# ########################## PLOTTING ############################


# # code to plot in a REPL window
# using Gtk
# c = Gtk.@GtkCanvas(400,300);
# w = Gtk.@GtkWindow(c,"data win");
# show(c);
# using Compose
# function sierpinski(n)
#     if n == 0
#         compose(context(), polygon([(1,1), (0,1), (1/2, 0)]))
#     else
#         t = sierpinski(n - 1)
#         compose(context(),
#                 (context(1/4,   0, 1/2, 1/2), t),
#                 (context(  0, 1/2, 1/2, 1/2), t),
#                 (context(1/2, 1/2, 1/2, 1/2), t))
#     end
# end
# co = sierpinski(5);
# Gtk.draw(c) do widget
#     Compose.draw(CAIROSURFACE(c.back),co)
# end
# Gtk.draw(c)
