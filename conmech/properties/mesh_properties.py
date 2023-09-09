from abc import ABC
from dataclasses import dataclass
from typing import List

import numpy as np
import meshio
from conmech.mesh.zoo.raw_mesh import RawMesh

@dataclass
class MeshDescription(ABC):
    initial_position: np.ndarray

    def build(self):
        raise NotImplementedError()


@dataclass
class ImportedMeshDescription(MeshDescription):
    path: str

    def build(self):
        mesh = meshio.read(self.path)
        return RawMesh(nodes=mesh.points, elements=mesh.cells_dict["triangle"])


@dataclass
class GeneratedMeshDescription(MeshDescription, ABC):
    max_element_perimeter: float


@dataclass
class HardcodedMeshDescription(MeshDescription, ABC):
    pass

@dataclass
class RectangleMeshDescription(GeneratedMeshDescription):
    scale: List[float]

    def build(self):
        from conmech.mesh.zoo.rectangle import Rectangle
        return Rectangle(self)


@dataclass
class CrossMeshDescription(GeneratedMeshDescription):
    scale: List[float]

    def build(self):
        from conmech.mesh.zoo.cross_for_tests import CrossMesh
        return CrossMesh(self)


@dataclass
class CubeMeshDescription(HardcodedMeshDescription):

    def build(self):
        from conmech.mesh.zoo.cube import Cube
        return Cube(self)


@dataclass
class BallMeshDescription(HardcodedMeshDescription):

    def build(self):
        from conmech.mesh.zoo.ball import Ball
        return Ball(self)


@dataclass
class Barboteu2008MeshDescription(GeneratedMeshDescription):
    
    def build(self):
        from conmech.mesh.zoo.barboteu_2008 import Barboteu2008
        return Barboteu2008(self)


@dataclass
class JOB2023MeshDescription(GeneratedMeshDescription):

    def build(self):
        from conmech.mesh.zoo.jurochbar_2023 import JOB2023
        return JOB2023(self)


@dataclass
class SOB2023MeshDescription(GeneratedMeshDescription):
    scale: List[float]

    def build(self):
        from conmech.mesh.zoo.sofochbar_2023 import SOB2023
        return SOB2023(self)


@dataclass
class CircleMeshDescription(GeneratedMeshDescription):
    radius: float

    def build(self):
        from conmech.mesh.zoo.pygmsh.dim_2.circle import Circle
        return Circle(self)


@dataclass
class PolygonMeshDescription(GeneratedMeshDescription):
    scale: List[float]

    def build(self):
        from conmech.mesh.zoo.pygmsh.dim_2.polygon import Polygon
        return Polygon(self)


@dataclass
class PgmshRectangleMeshDescription(GeneratedMeshDescription):
    scale: List[float]

    def build(self):
        from conmech.mesh.zoo.pygmsh.dim_2.rectangle import PgmshRectangle
        return PgmshRectangle(self)


@dataclass
class SplineMeshDescription(GeneratedMeshDescription):
    scale: List[float]

    def build(self):
        from conmech.mesh.zoo.pygmsh.dim_2.spline import Spline
        return Spline(self)


@dataclass
class Polygon3DMeshDescription(GeneratedMeshDescription):

    def build(self):
        from conmech.mesh.zoo.pygmsh.dim_3.polygon import Polygon3D
        return Polygon3D(self)


@dataclass
class TwistMeshDescription(GeneratedMeshDescription):

    def build(self):
        from conmech.mesh.zoo.pygmsh.dim_3.twist import Twist
        return Twist(self)
