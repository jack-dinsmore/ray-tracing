import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import convolve
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.colors import LinearSegmentedColormap

REMOVE_POLE = 1
CORONA_PROTON_GAMMA_MINUS_ONE = 0.1
CORONA_DENSITY_HEIGHT = 2.25 / CORONA_PROTON_GAMMA_MINUS_ONE
CORONA_ELECTRON_GAMMA = 1836.0 * CORONA_PROTON_GAMMA_MINUS_ONE
RESCALE_FOR_XRAY = CORONA_ELECTRON_GAMMA * CORONA_ELECTRON_GAMMA

X_RAY_WEIGHT = 10

bk = 8.61733326e-5

def to_xray(arr):
    res = arr[[1,2,0],:,:]
    res[1,:,:] /= 2
    return res

plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "cm"

figsize = (16, 9)
FALLOFF_SIG = 1000.0

def get_colormap(xray=False):
    def temp_to_color(temp):
        TEMP_TO_COLOR = [
            (255, 56, 0),
            (255, 71, 0),
            (255, 83, 0),
            (255, 93, 0),
            (255, 101, 0),
            (255, 109, 0),
            (255, 115, 0),
            (255, 121, 0),
            (255, 126, 0),
            (255, 131, 0),
            (255, 138, 18),
            (255, 142, 33),
            (255, 147, 44),
            (255, 152, 54),
            (255, 157, 63),
            (255, 161, 72),
            (255, 165, 79),
            (255, 169, 87),
            (255, 173, 94),
            (255, 177, 101),
            (255, 180, 107),
            (255, 184, 114),
            (255, 187, 120),
            (255, 190, 126),
            (255, 193, 132),
            (255, 196, 137),
            (255, 199, 143),
            (255, 201, 148),
            (255, 204, 153),
            (255, 206, 159),
            (255, 209, 163),
            (255, 211, 168),
            (255, 213, 173),
            (255, 215, 177),
            (255, 217, 182),
            (255, 219, 186),
            (255, 221, 190),
            (255, 223, 194),
            (255, 225, 198),
            (255, 227, 202),
            (255, 228, 206),
            (255, 230, 210),
            (255, 232, 213),
            (255, 233, 217),
            (255, 235, 220),
            (255, 236, 224),
            (255, 238, 227),
            (255, 239, 230),
            (255, 240, 233),
            (255, 242, 236),
            (255, 243, 239),
            (255, 244, 242),
            (255, 245, 245),
            (255, 246, 247),
            (255, 248, 251),
            (255, 249, 253),
            (254, 249, 255),
            (252, 247, 255),
            (249, 246, 255),
            (247, 245, 255),
            (245, 243, 255),
            (243, 242, 255),
            (240, 241, 255),
            (239, 240, 255),
            (237, 239, 255),
            (235, 238, 255),
            (233, 237, 255),
            (231, 236, 255),
            (230, 235, 255),
            (228, 234, 255),
            (227, 233, 255),
            (225, 232, 255),
            (224, 231, 255),
            (222, 230, 255),
            (221, 230, 255),
            (220, 229, 255),
            (218, 229, 255),
            (217, 227, 255),
            (216, 227, 255),
            (215, 226, 255),
            (214, 225, 255),
            (212, 225, 255),
            (211, 224, 255),
            (210, 223, 255),
            (209, 223, 255),
            (208, 222, 255),
            (207, 221, 255),
            (207, 221, 255),
            (206, 220, 255),
            (205, 220, 255),
            (207, 218, 255),
            (207, 218, 255),
            (206, 217, 255),
            (205, 217, 255),
            (204, 216, 255),
            (204, 216, 255),
            (203, 215, 255),
            (202, 215, 255),
            (202, 214, 255),
            (201, 214, 255),
            (200, 213, 255),
            (200, 213, 255),
            (199, 212, 255),
            (198, 212, 255),
            (198, 212, 255),
            (197, 211, 255),
            (197, 211, 255),
            (197, 210, 255),
            (196, 210, 255),
            (195, 210, 255),
            (195, 209, 255)
        ]

        float_index = temp / 100.0 - 10.0
        if float_index < 0.0:
            falloff = np.exp(-float_index * float_index / FALLOFF_SIG)
            return (
                TEMP_TO_COLOR[0][0] / 255.0 * falloff,
                TEMP_TO_COLOR[0][1] / 255.0 * falloff,
                TEMP_TO_COLOR[0][2] / 255.0 * falloff
            )
        elif float_index > len(TEMP_TO_COLOR) - 2.0:
            falloff = np.exp(-(float_index - len(TEMP_TO_COLOR) - 1.0)**2 / FALLOFF_SIG)
            return (
                TEMP_TO_COLOR[len(TEMP_TO_COLOR)- 1][0] / 255.0 * falloff,
                TEMP_TO_COLOR[len(TEMP_TO_COLOR)- 1][1] / 255.0 * falloff,
                TEMP_TO_COLOR[len(TEMP_TO_COLOR)- 1][2] / 255.0 * falloff
            )
        
        int_index = int(float_index)
        frac = float_index - int_index
        return [
            (TEMP_TO_COLOR[int_index + 1][0] * frac + TEMP_TO_COLOR[int_index][0] * (1.0 - frac)) / 255.0,
            (TEMP_TO_COLOR[int_index + 1][1] * frac + TEMP_TO_COLOR[int_index][1] * (1.0 - frac)) / 255.0,
            (TEMP_TO_COLOR[int_index + 1][2] * frac + TEMP_TO_COLOR[int_index][2] * (1.0 - frac)) / 255.0,
        ]

    cdict = {'red': [], 'green': [], 'blue': [] }
    min_temp = 1000
    max_temp = (111 + 10) * 100 / 2
    temps = np.arange(min_temp, max_temp,100)
    indices = np.linspace(0, 1, len(temps))
    for i, temp in zip(indices, temps):
        if xray:
            r, g, b = to_xray(np.array(temp_to_color(temp)).reshape(3,1,1))[:,0,0]
        else:
            r, g, b = temp_to_color(temp)
        if i != 0:
            cdict['red'][-1][-1] = r
            cdict['green'][-1][-1] = g
            cdict['blue'][-1][-1] = b
        cdict['red'].append([i, r, r])
        cdict['green'].append([i, g, g])
        cdict['blue'].append([i, b, b])
    blue_green1 = LinearSegmentedColormap('thermalmap', cdict)
    plt.figure()

    if xray:
        scale = RESCALE_FOR_XRAY/1000
    else:
        scale = 1

    c = plt.imshow([[min_temp * bk * scale, max_temp * bk * scale]], cmap=blue_green1)
    return c

def process(tag, max_flux=100, extra_name=""):
    with open(f"../data/{tag}-optical.npy", 'rb') as f:
        optical = np.load(f)

    with open(f"../data/{tag}-xray.npy", 'rb') as f:
        xray = np.load(f)

    colors = [optical, to_xray(xray)]

    for i in range(len(colors)):
        colors[i] = np.transpose(colors[i], axes=(1,2,0))
        colors[i] = np.flip(colors[i], axis=0)
        max_val = np.nanpercentile(colors[i], max_flux)
        for x,y in np.transpose(np.where(np.any(np.isnan(colors[i]), axis=-1))):
            if x < colors[i].shape[0] - 1 and y < colors[i].shape[1] - 1:
                colors[i][x,y,:] = np.nanmean([
                    colors[i][x,y+1,:],
                    colors[i][x,y-1,:],
                    colors[i][x+1,y,:],
                    colors[i][x-1,y,:],
                ],axis=0)
        colors[i] /= np.maximum(max_val, 1)
        if REMOVE_POLE is not None:
            middle = colors[i].shape[1]//2
            true_column = colors[i][:, middle - REMOVE_POLE - 1] / 2
            true_column += colors[i][:, middle + REMOVE_POLE + 1] / 2
            for j in range(middle - REMOVE_POLE, middle + REMOVE_POLE + 1):
                colors[i][:,j] = true_column

    fig, ax = plt.subplots(figsize=figsize, frameon=False)
    ax.set_axis_off()
    ax.imshow(colors[0], vmin=0, vmax=1)
    axins = inset_axes(ax, width="20%", height="1.3%", loc='upper left', borderpad=1)
    cbar = plt.colorbar(get_colormap(), cax=axins, orientation='horizontal')
    cbar_xticks = plt.getp(cbar.ax.axes, 'xticklabels')
    cbar.set_label("Optical Energy (eV)", color="white", usetex=True)
    plt.setp(cbar_xticks, color="white")
    ax.set_aspect("equal")
    fig.savefig(f"../data/{tag}-{extra_name}optical.png", bbox_inches='tight', pad_inches=0)
    fig.savefig(f"../data/{tag}-{extra_name}optical.pdf", bbox_inches='tight', pad_inches=0)

    if np.nanmax(colors[1][np.isfinite(colors[1])]) > 1:
        fig, ax = plt.subplots(figsize=figsize, frameon=False)
        ax.set_axis_off()
        ax.imshow(colors[1], vmin=0, vmax=1)
        axins = inset_axes(ax, width="20%", height="1.3%", loc='upper left', borderpad=1)
        cbar = plt.colorbar(get_colormap(xray=True), cax=axins, orientation='horizontal')
        cbar_xticks = plt.getp(cbar.ax.axes, 'xticklabels')
        cbar.set_label("X-ray Energy (keV)", color="white", usetex=True)
        plt.setp(cbar_xticks, color="white")
        ax.set_aspect("equal")
        fig.savefig(f"../data/{tag}-{extra_name}xray.png", bbox_inches='tight', pad_inches=0)
        fig.savefig(f"../data/{tag}-{extra_name}xray.pdf", bbox_inches='tight', pad_inches=0)

def together(tag, max_flux=100):
    with open(f"../data/{tag}-optical.npy", 'rb') as f:
        optical = np.load(f)

    with open(f"../data/{tag}-xray.npy", 'rb') as f:
        xray = np.load(f)

    colors = optical + X_RAY_WEIGHT * to_xray(xray)

    colors = np.transpose(colors, axes=(1,2,0))
    colors = np.flip(colors, axis=0)
    max_val = np.nanpercentile(colors, max_flux)
    for x,y in np.transpose(np.where(np.any(np.isnan(colors), axis=-1))):
        if x < colors.shape[0] - 1 and y < colors.shape[1] - 1:
            colors[x,y,:] = np.nanmean([
                colors[x,y+1,:],
                colors[x,y-1,:],
                colors[x+1,y,:],
                colors[x-1,y,:],
            ],axis=0)
    colors /= np.maximum(max_val, 1)

    if REMOVE_POLE is not None:
        middle = colors.shape[1]//2
        true_column = colors[:, middle - REMOVE_POLE - 1] / 2
        true_column += colors[:, middle + REMOVE_POLE + 1] / 2
        for j in range(middle - REMOVE_POLE, middle + REMOVE_POLE + 1):
                colors[:,j] = true_column

    fig, ax = plt.subplots(figsize=figsize, frameon=False)
    ax.set_axis_off()
    ax.imshow(colors, vmin=0, vmax=1)

    axins = inset_axes(ax, width="20%", height="1.3%", loc='upper left', borderpad=1)
    cbar = plt.colorbar(get_colormap(), cax=axins, orientation='horizontal')
    cbar_xticks = plt.getp(cbar.ax.axes, 'xticklabels')
    cbar.set_label("Optical Energy (eV)", color="white", usetex=True)
    plt.setp(cbar_xticks, color="white")

    axins = inset_axes(ax, width="20%", height="1.3%",
        bbox_to_anchor=(0.0,-0.008,1,1), bbox_transform=ax.transAxes)
    cbar = plt.colorbar(get_colormap(xray=True), cax=axins, orientation='horizontal')
    cbar_xticks = plt.getp(cbar.ax.axes, 'xticklabels')
    cbar.set_label("X-ray Energy (keV)", color="white", usetex=True)
    plt.setp(cbar_xticks, color="white")

    ax.set_aspect("equal")
    fig.savefig(f"../data/{tag}-tog.png", bbox_inches='tight', pad_inches=0)
    fig.savefig(f"../data/{tag}-tog.pdf", bbox_inches='tight', pad_inches=0)


if __name__ == "__main__":
    low_flux = 90
    high_flux = 99.5
    process("flat", max_flux=low_flux)
    process("minkowski", max_flux=low_flux)
    process("thick", max_flux=low_flux)
    process("thick", max_flux=high_flux, extra_name='bright-')
    process("thin", max_flux=low_flux)
    process("thin", max_flux=high_flux, extra_name='bright-')
    process("schwarzschild", max_flux=low_flux)
    process("schwarzschild", max_flux=high_flux, extra_name='bright-')
    process("kerr", max_flux=low_flux)
    process("kerr", max_flux=high_flux, extra_name='bright-')

    together("schwarzschild", max_flux=90)
    together("kerr", max_flux=98)
