#!/usr/bin/env ovitos
import numpy as np
import sys
import ovito
import re
import argparse

pipeline=ovito.pipeline.Pipeline()
pipeline.source=ovito.pipeline.FileSource()
cna_mod=ovito.modifiers.CommonNeighborAnalysisModifier(mode=ovito.modifiers.CommonNeighborAnalysisModifier.Mode.AdaptiveCutoff)
cna_mod.structures[ovito.modifiers.CommonNeighborAnalysisModifier.Type.FCC].color =(223/255,223/255,223/255)
cna_mod.structures[ovito.modifiers.CommonNeighborAnalysisModifier.Type.BCC].color =(102/255,102/255,255/255)
cna_mod.structures[ovito.modifiers.CommonNeighborAnalysisModifier.Type.HCP].color =(255/255,88/255,88/255)
cna_mod.structures[ovito.modifiers.CommonNeighborAnalysisModifier.Type.OTHER].color =(78/255,246/255,255/255)
pipeline.modifiers.append(cna_mod)

# exp_select_mod=ovito.modifiers.ExpressionSelectionModifier(expression = 'v_group_label==6 || v_group_label==7')
# asg_color_mod=ovito.modifiers.AssignColorModifier(color=(1,0.333,0.498))
# pipeline.modifiers.append(exp_select_mod)
# pipeline.modifiers.append(asg_color_mod)
# exp_select_mod=ovito.modifiers.ExpressionSelectionModifier(expression = 'v_group_label==4 || v_group_label==5')
# asg_color_mod=ovito.modifiers.AssignColorModifier(color=(1,1,0.498))
# pipeline.modifiers.append(exp_select_mod)
# pipeline.modifiers.append(asg_color_mod)

pipeline.add_to_scene()
vp=ovito.vis.Viewport()
vp.type=ovito.vis.Viewport.Type.Ortho
vp.camera_dir=(0,0,-1) # vp.camera_up=(0,1,0) # not recognized

full={'fov':320,'size':(1000,1000),'radius':0.8}
# enlarge={'fov':60,'size':(400,600),'radius':0.8}
enlarge={'fov':50,'size':(160,240),'radius':0.8}
tip={'fov':8,'size':(100,100),'radius':0.6}
# custom={'fov':6,'size':(400,60)}
custom={'fov':8,'size':(100,75),'radius':0.6}
ir={'fov':40,'size':(160,200),'radius':0.8}
mode_opt={'full':full,'enlarge':enlarge,'tip':tip,'custom':custom,'ir':ir}

def crack_pos(fname='./data/crack_pos.mod'):
  regex=re.compile('variable\s+(Xc|Yc)\s+equal\s+(\S+)')
  with open(fname,'r') as f:
    for line in f:
      mo=regex.search(line)
      if mo:
        if mo[1]=='Xc':
          Xc=float(mo[2])
        elif mo[1]=='Yc':
          Yc=float(mo[2])
  return (Xc, Yc)
def read_k_from_fname(fname):
  regex=re.compile('dump.k_load_(-?\d+\.\d+)_(-?\d+\.\d+)_(-?\d+\.\d+)')
  mo=regex.search(fname)
  if mo:
    kii=float(mo[1]); ki=float(mo[2]); kiii=float(mo[3])
  else:
    raise Exception("Wrong fname format")
  return (kii,ki,kiii)
def text_overlay(text,x=0,y=0,fontsize=0.03):
    overlay = ovito.vis.TextLabelOverlay(
        text = text,
        # alignment = ovito.qt_compat.QtCore.Qt.AlignmentFlag.AlignLeft | ovito.qt_compat.QtCore.Qt.AlignmentFlag.AlignTop,
        offset_x = x,
        offset_y = y,
        font_size = fontsize,
        text_color = (0,0,0))
    return overlay
def line_overlay(x1,y1,x2,y2):
  def pyol(args: ovito.vis.PythonViewportOverlay.Arguments):
    x1p,y1p=args.project_point((x1,y1,0))
    x2p,y2p=args.project_point((x2,y2,0))
    args.painter.drawLine(x1p,y1p,x2p,y2p)
  overlay=ovito.vis.PythonViewportOverlay(function=pyol)
  return overlay
def apply_overlay():
  # kii,ki,kiii=read_k_from_fname(fname)
  # vp.overlays.append(text_overlay("{ki:.2f}".format(ki=ki),0.7,0))
  # add cross mark at crack tip
  length=0.5; Xc, Yc=crack_pos();
  vp.overlays.append(line_overlay(Xc,Yc-length,Xc,Yc+length))
  vp.overlays.append(line_overlay(Xc-length,Yc,Xc+length,Yc))
def render(fname,ofile,mode):
  pipeline.source.load(fname)
  data=pipeline.compute()
  data.cell.vis.enabled=False
  # data.particles.vis.radius=0.8
  data.particles.vis.radius=mode_opt[mode]['radius']
  Xc, Yc=crack_pos(); vp.camera_pos=(Xc,Yc,0)
  vp.overlays=[]
  apply_overlay()
  vp.fov=mode_opt[mode]['fov']
  ofname='./figure/'+fname+ofile+'.png'
  vp.render_image(size=mode_opt[mode]['size'],filename=ofname,renderer=ovito.vis.OpenGLRenderer())

def main():
   # ovitos ovito_visualization.py -f fname1 fname2 fname3 ...
  parser=argparse.ArgumentParser(
    prog="python3 ovito_visualization.py",
    description="Render image for an atomistic model",
    epilog="Written by LI Songwei")
  parser.add_argument('-f','--fnames',type=str,required=True,nargs='*',help="Dump files to be rendered.")
  parser.add_argument('-o','--ofile',type=str,default='',help='Output file name.')
  parser.add_argument('-m','--mode',type=str,default='enlarge',choices=['full','enlarge','tip','custom','ir'],help='Region to be captured.')
  args=parser.parse_args()
  for f in args.fnames:
    render(f,args.ofile,args.mode)

if __name__ == "__main__":
  main()
