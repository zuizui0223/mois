#!/usr/bin/env python3
"""One-pass floral-region dominant colour in device-sRGB-derived CIELAB.

The 500px YOLO FLOWER bounding box is a pragmatic proxy for visible petal/
corolla colour across morphological groups, not an annotated petal mask.
Explicitly retain green and brown-like colours: neither is reclassified as
"white" or dropped as background automatically.

Only samples bounded 80 x 80 RGB pixels (already downloaded in Stage 1).
No EfficientSAM, semantic organ classifier, CIELAB-from-nine-palette guessing,
image re-download, raw image persistence or calibrated reflectance claims.

Protocol v1: D65 sRGB conversion, 12.5-L* / 20-a* / 20-b* modal bins,
median modal-bin Lab, robust continuous L*/C* medians, and chroma-weighted
circular hue. Low-chroma hue is UNDEFINED. Fixed semantic-looking fractions
(white-like, green-like, brown-like) are candidate-pixel descriptors and
NOT botanical ground-truth classification. Do not exclude any bins by colour.
"""
from __future__ import annotations

import math
import numpy as np

PROTOCOL="mois_flower_box_srgb_cielab_mode_v1"
SIZE=(80,80)
DESCRIPTOR_COLUMNS=(
    "floral_box_colour_protocol",
    "floral_box_main_L_star","floral_box_main_a_star",
    "floral_box_main_b_star","floral_box_main_C_star",
    "floral_box_main_hue_degrees",
    "floral_box_main_mode_fraction","floral_box_second_mode_fraction",
    "floral_box_main_rgb_hex",
    "floral_box_L50","floral_box_L10","floral_box_L90",
    "floral_box_C50","floral_box_C10","floral_box_C90",
    "floral_box_hue_circular_mean_degrees",
    "floral_box_fraction_low_chroma",
    "floral_box_fraction_white_like",
    "floral_box_fraction_green_like",
    "floral_box_fraction_brown_like",
    "floral_box_central_vs_border_lab_delta_e",
    "floral_box_may_include_nonfloral_pixels",
)
def rgb_to_lab(rgb):
    """Convert RGB uint8 or [0,1] sRGB to D65 CIE L*a*b* using numpy."""
    arr=np.asarray(rgb)
    if arr.ndim<1 or arr.shape[-1]!=3:
        raise ValueError("RGB array must have last dimension 3")
    if np.issubdtype(arr.dtype,np.integer):
        x=np.asarray(arr,dtype=np.float64)/255.0
    else:
        x=np.asarray(arr,dtype=np.float64)
    if not np.isfinite(x).all() or np.any((x<0)|(x>1)):
        raise ValueError("sRGB values must be finite in [0,1]")
    linear=np.where(x<=.04045,x/12.92,((x+.055)/1.055)**2.4)
    xyz=linear@np.array([
        [.4124564,.3575761,.1804375],
        [.2126729,.7151522,.0721750],
        [.0193339,.1191920,.9503041]
    ]).T
    norm=xyz/np.array([.95047,1.0,1.08883])
    delta=6.0/29.0
    f=np.where(norm>delta**3,np.cbrt(np.maximum(norm,0)),
               norm/(3*delta**2)+4.0/29.0)
    return np.stack([
        116*f[...,1]-16,
        500*(f[...,0]-f[...,1]),
        200*(f[...,1]-f[...,2]),
    ],axis=-1)

def circular_mean_degrees(angles,weights):
    angles=np.asarray(angles,dtype=float)
    w=np.asarray(weights,dtype=float)
    if not len(angles) or w.sum()<=0:
        return None
    theta=np.deg2rad(angles)
    x=np.sum(np.cos(theta)*w)
    y=np.sum(np.sin(theta)*w)
    if np.hypot(x,y)<=1e-8*w.sum():
        return None
    return float(np.degrees(np.arctan2(y,x))%360.0)

def mode_of_lab(lab,rgb):
    """O(N) histogram of quantized CIELAB; one representative dominant mode."""
    flat=lab.reshape((-1,3))
    rgbflat=rgb.reshape((-1,3))
    bins=np.column_stack((
        np.floor(np.clip(flat[:,0],0,99.999999)/12.5).astype(np.int32),
        np.floor((flat[:,1]+128)/20).astype(np.int32),
        np.floor((flat[:,2]+128)/20).astype(np.int32),
    ))
    indexes,rev,count=np.unique(bins,axis=0,return_inverse=True,return_counts=True)
    # Ties resolved consistently by lexicographically ordered Lab bin.
    ranking=np.argsort(-count,kind="stable")
    first=(rev==ranking[0])
    primary=np.median(flat[first],axis=0)
    representative_rgb=np.median(rgbflat[first],axis=0).astype(np.uint8)
    second_ratio=float(count[ranking[1]]/len(flat)) if len(ranking)>1 else 0.
    return primary,float(count[ranking[0]]/len(flat)),second_ratio,(
        "#"+"".join(f"{int(v):02X}" for v in representative_rgb))

def _format(value,digits=5):
    if value is None:
        return ""
    v=float(value)
    return f"{v:.{digits}f}" if math.isfinite(v) else ""

def flower_box_lab_descriptor(image,xyxy):
    """Measured pixel colour features inside a YOLO flower box.

    The 'major petal colour' is ONLY an operational estimate of the dominant
    colour in this detected floral ROI, not separated anatomical petal truth.
    """
    from PIL import Image
    if not isinstance(image,Image.Image):
        raise ValueError("Expected decoded Pillow RGB image")
    if len(xyxy)!=4 or not all(math.isfinite(float(v)) for v in xyxy):
        raise ValueError("Invalid YOLO flower box geometry")
    a,b,c,d=map(float,xyxy)
    l=max(0,min(image.width-1,int(math.floor(a))))
    t=max(0,min(image.height-1,int(math.floor(b))))
    r=max(l+1,min(image.width,int(math.ceil(c))))
    bot=max(t+1,min(image.height,int(math.ceil(d))))
    rgb=np.asarray(image.convert("RGB").crop((l,t,r,bot)).resize(
        SIZE,Image.Resampling.BILINEAR),dtype=np.uint8)
    lab=rgb_to_lab(rgb)
    L=lab[:,:,0]
    aa=lab[:,:,1]
    bb=lab[:,:,2]
    C=np.hypot(aa,bb)
    hue=np.degrees(np.arctan2(bb,aa))%360.
    dom,first,second,hexcode=mode_of_lab(lab,rgb)
    main_C=float(np.hypot(dom[1],dom[2]))
    main_hue=float(np.degrees(np.arctan2(dom[2],dom[1]))%360) if main_C>=12 else None
    hue_valid=C>=12
    hmean=circular_mean_degrees(hue[hue_valid],C[hue_valid])
    # These are reported as explicit thumbnail BOX-PIXEL diagnostics, not
    # petal colour frequencies or independent flower biological categories.
    white=(L>=80)&(C<=15)
    green=(C>=20)&(hue>=105)&(hue<=190)
    brown=(C>=18)&(L<55)&(hue>=20)&(hue<=100)
    # Central-vs-border deltaE: uncertainty indicator, never remove greens.
    h,w=L.shape
    yy,xx=np.indices((h,w))
    middle=((yy>=h//4)&(yy<3*h//4)&
            (xx>=w//4)&(xx<3*w//4))
    center=np.median(lab[middle],axis=0)
    border=np.median(lab[~middle],axis=0)
    delta=float(np.linalg.norm(center-border))
    metrics={
        "floral_box_colour_protocol":PROTOCOL,
        "floral_box_main_L_star":_format(dom[0]),
        "floral_box_main_a_star":_format(dom[1]),
        "floral_box_main_b_star":_format(dom[2]),
        "floral_box_main_C_star":_format(main_C),
        "floral_box_main_hue_degrees":_format(main_hue),
        "floral_box_main_mode_fraction":_format(first,6),
        "floral_box_second_mode_fraction":_format(second,6),
        "floral_box_main_rgb_hex":hexcode,
        "floral_box_L50":_format(np.median(L)),
        "floral_box_L10":_format(np.quantile(L,.1)),
        "floral_box_L90":_format(np.quantile(L,.9)),
        "floral_box_C50":_format(np.median(C)),
        "floral_box_C10":_format(np.quantile(C,.1)),
        "floral_box_C90":_format(np.quantile(C,.9)),
        "floral_box_hue_circular_mean_degrees":_format(hmean),
        "floral_box_fraction_low_chroma":_format(np.mean(C<12),6),
        "floral_box_fraction_white_like":_format(np.mean(white),6),
        "floral_box_fraction_green_like":_format(np.mean(green),6),
        "floral_box_fraction_brown_like":_format(np.mean(brown),6),
        "floral_box_central_vs_border_lab_delta_e":_format(delta),
        "floral_box_may_include_nonfloral_pixels":"true",
    }
    if set(metrics)!=set(DESCRIPTOR_COLUMNS):
        raise RuntimeError("Unversioned CIELAB descriptor schema drift")
    return metrics
