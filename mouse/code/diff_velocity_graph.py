import math

def percentile(data, perc: int):
    size = len(data)
    return sorted(data)[int(math.ceil((size * perc) / 100)) - 1]

adatas = adatas['ko']
def diff_velocity_graph(adatas, basis, label, components=None, vkey='velocity',
legend_loc='none', legend_fontsize=None, legend_fontweight=None,
refkey='ko', altkey='wt'):
    keys = [refkey, altkey]
    keys.sort()
    adatakeys = list(adatas.keys())
    adatakeys.sort()
    if adatakeys != keys:
        raise ValueError('refkey and altkey dont match adatas keys')
    
    X_grids = {}
    V_grids = {}
    for adata_id in list(adatas.keys()):
        adata = adatas[adata_id]
        comps, obsm = scv.pl.utils.get_components(components, basis), adata.obsm
        X_emb = np.array(obsm[f"X_{basis}"][:, comps])
        V_emb = np.array(obsm[f"{vkey}_{basis}"][:, comps])
        n_obs, n_dim = X_emb.shape
        X_grid, V_grid = compute_velocity_on_grid(
            X_emb=X_emb,
            V_emb=V_emb,
            density=None,
            autoscale=True,
            smooth=None,
            n_neighbors=None,
            min_mass=None,
        )
        X_grids[adata_id] = X_grid
        V_grids[adata_id] = V_grid
    
    # Match the basis_reduction maps between KO and WT conditions to allow for
    # velocity comparison of equivalent locations
    maxdist = 0.05
    match_arr = np.empty((0,2))
    for row_ko in X_grids[refkey]:
        row_ko_match=[]
        for row_wt in X_grids[altkey]:
            # absolute distance between ref_i and alt_J
            row_ko_match.append(np.absolute(np.array(row_ko) - np.array(row_wt)).sum())
        # For each row in REF, find the ALT row with the minimum distance
        minIndex, minValue = min(enumerate(row_ko_match), key=lambda v: v[1])
        minRow = [int(minIndex), round(minValue,2)]
        # match_arr is a 2-column matrix:
        #  col1 = for i = 1..I; row index from X_grids[altkey] that matches the row for X_grids[refkey][i,]
        #  col2 = for i = 1..I; the absolute difference between X_grids[altkey] row and X_grids[refkey][i,]
        match_arr = np.vstack([match_arr, minRow])
    
    # Convert X_grids[altkey] index to an int
    idx = list(match_arr[:,0])
    idx_int = list()
    for i in idx:
        idx_int.append(int(i))
    
    # Calculate the velocity difference between the best-index matched spots
    V_grid_unity = V_grids[refkey] - V_grids[altkey][idx_int,:]
    X_grid_unity = X_grids[refkey] # X_grids[altkey][idx_int,:]
    v_arr_sums = list()
    for v_arr in V_grid_unity:
        v_arr_sums.append(np.absolute(v_arr).sum())
    
    high_diff_idx =  np.nonzero([float(i) > percentile(list(v_arr_sums), 90) for i in v_arr_sums])
    high_diff_idx = list(high_diff_idx[0])
    
    
    
    #
    colors=label # 'manual_anno'
    color=colors
    multikey = (
        colors
        if len(colors) > 1
        else layers
        if len(layers) > 1
        else vkeys
        if len(vkeys) > 1
        else None
    )
    
    ax , show , figsize , dpi , arrow_size , arrow_length , arrow_color, scale = None, None, None, None, None, None, None, None
    ax, show = scv.pl.utils.get_ax(ax, show, figsize, dpi)
    hl, hw, hal = scv.pl.utils.default_arrow(arrow_size)
    if arrow_length is not None:
        scale = 1 / arrow_length
    
    if scale is None:
        scale = 1
    
    if arrow_color is None:
        arrow_color = "grey"
    
    quiver_kwargs = {"angles": "xy", "scale_units": "xy", "edgecolors": "k"}
    quiver_kwargs.update({"scale": scale, "width": 0.001, "headlength": hl / 2})
    quiver_kwargs.update({"headwidth": hw / 2, "headaxislength": hal / 2})
    quiver_kwargs.update({"color": arrow_color, "linewidth": 0.2, "zorder": 3})
    size = 4 * scv.pl.utils.default_size(adata)
    title = refkey.upper() + ' [red] - ' + altkey.upper() + ' [blue] difference'
    
    layer=None
    scatter_kwargs = {
        "basis": basis,
        "perc": None,
        "use_raw": None,
        "sort_order": True,
        "alpha": 0.2,
        "components": None,
        "projection": "2d",
        "legend_loc": legend_loc,
        "groups": None,
        "legend_fontsize": legend_fontsize,
        "legend_fontweight": legend_fontweight,
        "palette": None,
        "color_map": None,
        "frameon": None,
        "xlabel": None,
        "ylabel": None,
        "colorbar": True,
        "dpi": None,
        "fontsize": None,
        "show": False,
        "save": False,
    }
    
    
    #quiver_kwargs.update({"color": 'white', "headwidth" : hw*4, "headaxislength" : hal*4, "headlength" : hl*4})
    #ax.quiver(
    #    X_grid_unity_high[:, 0], X_grid_unity_high[:, 1], V_grid_unity_high[:, 0], V_grid_unity_high[:, 1], **quiver_kwargs
    #)
    quiver_kwargs.update({"color": 'grey'})
    ax.quiver(
        X_grids[refkey][:, 0], X_grids[refkey][:, 1], V_grids[refkey][:, 0], V_grids[refkey][:, 1], **quiver_kwargs
    )
    quiver_kwargs.update({"color": 'red', "headwidth" : hw*10, "headaxislength" : hal*10, "headlength" : hl*10})
    ax.quiver(
        X_grids[refkey][high_diff_idx, 0], X_grids[refkey][high_diff_idx, 1],
        V_grids[refkey][high_diff_idx, 0], V_grids[refkey][high_diff_idx, 1], **quiver_kwargs
    )
    quiver_kwargs.update({"color": 'blue', "headwidth" : hw*10, "headaxislength" : hal*10, "headlength" : hl*10})
    ax.quiver(
        X_grids[altkey][idx_int,:][high_diff_idx, 0], X_grids[altkey][idx_int,:][high_diff_idx, 1],
        V_grids[altkey][idx_int,:][high_diff_idx, 0], V_grids[altkey][idx_int,:][high_diff_idx, 1], **quiver_kwargs
    )
    
    
    ax = scv.pl.scatter(
        adata,
        layer=layer,
        color=color,
        size=size,
        title=title,
        ax=ax,
        zorder=0,
        **scatter_kwargs,
    )
    return(ax)
