using GLMakie
funcs = [1,2]
labels = ["theme_dark", "theme_ggplot2"]
theme_indx = Observable(1)

fig, ax, obj= scatter(1:10)
menu = Menu(fig[1, 1], options = zip(labels, funcs))
fig[1, 1, Left()] = vgrid!(
    Label(fig, "theme", width = 100),
    menu)
on(menu.selection) do s
    @show s
    theme_indx[] = s
    if s==1
        @show s
        set_theme!()
        set_theme!(theme_dark())
        notify(theme_indx)
    elseif s==2
        @show s
        set_theme!()
        set_theme!(theme_ggplot2())
        notify(theme_indx)
    end
end
fig
