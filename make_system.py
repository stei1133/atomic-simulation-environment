from ase import Atoms
from ase.build import add_adsorbate, add_vacuum, graphene_nanoribbon, mx2, sort, make_supercell, root_surface
from ase.io import read, write
from ase.spacegroup import crystal
from ase.visualize import view
import math
def make_system(name,h,a,n,m):
    print('running...')
    """Function that takes in the height, h, the spacing of the adsorbate, a,
    the number of positive pole down molecules, n, and the number of negative
    pole down molecules, m."""
    """name specifies the type of base layer: graphene or mos2"""
    d=0.92791085 #this value is from the relaxed HF molecule in a vacuum
    mol=Atoms('HF',positions=[(0,0,0),(0,0,d)])
    """mol is the dipole molecule. d indicates the bond length"""
    if name=='graphene':
        num=a*(n+m)/4.26
        """number of unit cells. 4.26 represents the cell length in angstroms"""
        global gnr
        gnr=graphene_nanoribbon(1,int(math.ceil(num)), type='armchair',vacuum=17.5,
                        saturated=False, sheet=True)
        """creates a graphene sheet with 10 A of vacuum on the top and bottom"""
        size=gnr.get_cell()
        cellwidth=size[0][0]
        cellheight=size[1][1]
        celllength=size[2][2]
        spacing=celllength/(n+m)
        """gets the cell dimensions"""
        mol.rotate(270,(1,0,0))
        for i in range(n):
            add_adsorbate(gnr, mol, height=(-i)*spacing,position=
                          (cellwidth/2,cellheight
                            /2+h+d),mol_index=1)
        """adds n HF molecules at a separation h from the graphene"""
        mol.rotate(180,(1,0,0))
        for j in range(m):
            add_adsorbate(gnr, mol, height=-(j+n)*spacing,position=
                          (cellwidth/2,
                         cellheight/2+h),mol_index=1)
            """adds m HF molecules in the opposite orientation"""
        write('graph'+str(n)+'up'+str(m)+'down_center.cif',gnr)
        write('graph'+str(n)+'up'+str(m)+'down_center.pwi',gnr)
    
    elif name=='BN':
        num=a*(n+m)/3
        global BN
        BN=Atoms('BN',positions=[(0,1.405056151,12.5),(1.2561747111,0.7252528076,12.5)])
        dim=[(2.5123494221,0,0),(-1.256174711,2.1757584227,0),(0,0,25)]
        BN.set_cell(dim)
        BN=BN.repeat((int(math.ceil(num)),2,1))
        dimensions=BN.get_cell()
        BN.set_cell([dimensions[0][0], (dimensions[1][1]),
                    dimensions[2][2]])
        cellwidth=dimensions[0][0]
        spacing=cellwidth/(n+m)
        for i in range(n):

            add_adsorbate(BN,mol,height=h,position=(i*spacing,
                                       dimensions[1][1]/2),mol_index=0)
        mol.rotate(180,(1,0,0))
        for j in range(m):
            add_adsorbate(BN,mol,height=h+d,position=((j+n)*spacing,
                                       dimensions[1][1]/2),mol_index=0)
        BN=sort(BN)
        write('BN.cif',BN)
        write(str(n)+'up'+str(m)+'down.pwi',BN)
    elif name=='mos2':
        """#this code runs for the MoS2 system
        num=a*(n+m)/3.17
        global mos2
        mos2=mx2(formula='MoS2', kind='2H', a=3.18, thickness=1, size=(int(math.ceil(num)),1, 1),
             vacuum=10)
        #mx2 makes an MoS2 system
        dim=mos2.get_cell()
        cellwidth=dim[0][0]
        cellheight=dim[1][1]
        celllength=dim[2][2]
        spacing=cellwidth/(n+m)
        #mos2.set_cell([(dim[0][0],0,0),(0,dim[1][1],0),(0,0,dim[2][2])],scale_Atoms=True)
        #truncates the unit cell to be cubic
        for i in range(n):
            add_adsorbate(mos2, mol, height=cellwidth+h,position=((-i)*spacing,cellheight/2
                                                                  ),mol_index=1)
            adds adsorbates to the MoS2
        mol.rotate(180,(1,0,0))
        for j in range(m):
            add_adsorbate(mos2, mol, height=cellwidth+h-d,position=(-(j+n)*spacing,
                                            cellheight/2),mol_index=1)
        write('mos2.cif',mos2)
        #write('mos2.pwi',mos2)"""
    print('complete')
    
make_system('graphene',5,3,2,2)
