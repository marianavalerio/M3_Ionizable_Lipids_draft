## Change the file names and run as specified in the line below:
in_file = 'vmdscene.stl'
out_file = 'out_test.tsi'
voxel_size = 0.1
# blender --background --python remesh.py

import bpy

### Remeshing the marching cubes mesh so it is watertight.

# Remove the initial cube
try:
    objs = bpy.data.objects
    objs.remove(objs["Cube"], do_unlink=True)
except:
    pass

# Import the mesh
bpy.ops.import_mesh.stl(filepath=in_file)

# Set the properties in the context
bpy.context.object.data.remesh_voxel_size = 0.1
bpy.context.object.data.use_remesh_fix_poles = True

# Remesh using default values
bpy.ops.object.voxel_remesh()

# Convert the quads to triangles
bpy.ops.object.editmode_toggle()
bpy.ops.mesh.select_all(action='SELECT')
bpy.ops.mesh.quads_convert_to_tris(quad_method='FIXED', ngon_method='CLIP')
bpy.ops.object.editmode_toggle()

### TS2CG TSI COVERTING PART
ob = bpy.context.active_object
bpy.ops.object.origin_set(type="GEOMETRY_ORIGIN")

box = []
for i in ob.dimensions:
    box.append(i*2)

f = open(out_file, 'w')
f.write('version 1.1\n')
f.write('box\t{}\t{}\t{}\n'.format(box[0], box[1], box[2]))
f.write('vertex\t{}\n'.format(len(ob.data.vertices)))

for v in ob.data.vertices:
    f.write('{}\t{:.6f}\t{:.6f}\t{:.6f}\n'.format(v.index, v.co.x+box[0]/2, v.co.y+box[1]/2, v.co.z+box[2]/2))

f.write('triangle\t{}\n'.format(len(ob.data.polygons)))
for t in ob.data.polygons:
    f.write('{}\t'.format(t.index))
    for v in t.vertices:
        f.write('{}\t'.format(v))
    f.write('\n')

f.write('inclusion\t0\n')
