# Constants for VGrid ID structure
HIERARCHY_LEVEL_BITS = 5  # Set to 5 bits (32 levels)
HIERARCHY_LEVEL_MASK = 2**HIERARCHY_LEVEL_BITS - 1
MAX_HIERARCHY_LEVEL = HIERARCHY_LEVEL_MASK

TILE_INDEX_BITS = 21  # Kept at 21 bits
TILE_INDEX_MASK = 2**TILE_INDEX_BITS - 1
MAX_TILE_INDEX = TILE_INDEX_MASK

OBJECT_INDEX_BITS = 21
OBJECT_INDEX_MASK = 2**OBJECT_INDEX_BITS - 1
MAX_OBJECT_INDEX = OBJECT_INDEX_MASK


class VGRID:
    """A custom grid class that allows defining cell sizes in degrees and resolution as integer."""

    def __init__(self, cell_size: float, aperture: int = 4):
        """
        Initialize a VGRID instance.

        Args:
            cell_size: Size of each cell in degrees
            aperture: Ratio to the next resolution (4 for quadtree-like, 9 for nonagon-like)
        """
        if cell_size <= 0:
            raise ValueError("cell_size must be positive")
        if 360 % cell_size != 0:
            raise ValueError("cell_size must divide 360 evenly")
        if 180 % cell_size != 0:
            raise ValueError("cell_size must divide 180 evenly")

        self.cell_size = cell_size
        self.aperture = aperture
        self.total_columns = int(360 / cell_size)
        self.total_rows = int(180 / cell_size)

        # Calculate resolution based on cell_size
        # Resolution 0: 1° cells (360 columns)
        # Resolution 1: 0.25° cells (1440 columns) for aperture=4
        # Resolution 2: 0.0625° cells (5760 columns) for aperture=4
        # etc.
        base_cell_size = 1.0  # 1° at resolution 0
        if cell_size == base_cell_size:
            self.resolution = 0
        else:
            # Calculate resolution: log_aperture(cell_size / base_cell_size)
            # For larger cell sizes, we need negative resolutions
            # For smaller cell sizes, we need positive resolutions
            import math

            if cell_size > base_cell_size:
                # Larger cells = lower resolution (negative)
                ratio = cell_size / base_cell_size
                self.resolution = -int(round(math.log(ratio, aperture)))
            else:
                # Smaller cells = higher resolution (positive)
                ratio = base_cell_size / cell_size
                self.resolution = int(round(math.log(ratio, aperture)))

        # Validate calculated resolution
        if self.resolution < -MAX_HIERARCHY_LEVEL:
            raise ValueError(
                f"Calculated resolution {self.resolution} is too negative (minimum: {-MAX_HIERARCHY_LEVEL})"
            )
        if self.resolution > MAX_HIERARCHY_LEVEL:
            raise ValueError(
                f"Calculated resolution {self.resolution} exceeds maximum {MAX_HIERARCHY_LEVEL}"
            )

    @classmethod
    def to_vgrid_id(cls, resolution: int, tile_index: int, object_index: int) -> int:
        """Create 64-bit representation of `VGridId`."""
        # Map negative resolutions to positive values in the upper half of the range
        if resolution < 0:
            stored_resolution = MAX_HIERARCHY_LEVEL + 1 + resolution
        else:
            stored_resolution = resolution

        assert 0 <= stored_resolution <= MAX_HIERARCHY_LEVEL
        assert 0 <= tile_index <= MAX_TILE_INDEX
        assert 0 <= object_index <= MAX_OBJECT_INDEX

        x = (
            # 5 bits
            stored_resolution
            # 21 bits, offset by previous field's width
            | (tile_index << HIERARCHY_LEVEL_BITS)
            # 21 bits, offset by the combined width of previous fields
            | (object_index << (HIERARCHY_LEVEL_BITS + TILE_INDEX_BITS))
        )
        return x

    @classmethod
    def get_resolution(cls, vgrid_id: int) -> int:
        """Resolution from 64-bit representation of `VGridId`."""
        stored_resolution = vgrid_id & HIERARCHY_LEVEL_MASK
        # Map back from stored resolution to actual resolution
        if stored_resolution > MAX_HIERARCHY_LEVEL // 2:
            # This was a negative resolution
            return -(MAX_HIERARCHY_LEVEL + 1 - stored_resolution)
        else:
            # This was a positive resolution
            return stored_resolution

    @classmethod
    def get_tile_index(cls, vgrid_id: int) -> int:
        """Tile index from 64-bit representation of `VGridId`."""
        offset = HIERARCHY_LEVEL_BITS
        return (vgrid_id >> offset) & TILE_INDEX_MASK

    @classmethod
    def get_object_index(cls, vgrid_id: int) -> int:
        """Object (node or edge) index from 64-bit representation of `VGridId`."""
        offset = HIERARCHY_LEVEL_BITS + TILE_INDEX_BITS
        return (vgrid_id >> offset) & OBJECT_INDEX_MASK

    def get_latlon(self, vgrid_id: int) -> tuple[float, float]:
        """Get latitude and longitude of tile's bottom-left corner."""
        resolution = self.get_resolution(vgrid_id)
        if resolution != self.resolution:
            raise ValueError(
                f"VGrid ID resolution {resolution} doesn't match VGRID instance resolution {self.resolution}"
            )

        tile_index = self.get_tile_index(vgrid_id)
        lat = (tile_index // self.total_columns) * self.cell_size - 90
        lon = (tile_index % self.total_columns) * self.cell_size - 180

        return lat, lon

    def get_tile_index_from_latlon(self, lat: float, lon: float) -> int:
        """Get tile index from latitude and longitude coordinates."""
        assert -90 <= lat <= 90, "Latitude must be between -90 and 90"
        assert -180 <= lon <= 180, "Longitude must be between -180 and 180"

        num_columns = int((lon + 180) / self.cell_size)
        num_rows = int((lat + 90) / self.cell_size)

        return num_rows * self.total_columns + num_columns

    def create_vgrid_id(self, lat: float, lon: float, object_index: int = 0) -> int:
        """Create VGrid ID from latitude, longitude and optional object index."""
        tile_index = self.get_tile_index_from_latlon(lat, lon)
        return self.to_vgrid_id(self.resolution, tile_index, object_index)

    def tiles_for_bounding_box(
        self,
        left: float,
        bottom: float,
        right: float,
        top: float,
        resolution: int = None,
    ) -> list[tuple[int, int]]:
        """Return a list of tiles that intersect the bounding box.

        Args:
            left: Left longitude boundary
            bottom: Bottom latitude boundary
            right: Right longitude boundary
            top: Top latitude boundary
            resolution: Resolution level to generate tiles at. If None, uses VGRID instance resolution.

        Returns:
            List of tuples (resolution, tile_index)
        """
        assert -90 <= bottom <= 90, "Bottom latitude must be between -90 and 90"
        assert -90 <= top <= 90, "Top latitude must be between -90 and 90"
        assert bottom <= top, (
            "Bottom latitude must be less than or equal to top latitude"
        )
        assert -180 <= left <= 180, "Left longitude must be between -180 and 180"
        assert -180 <= right <= 180, "Right longitude must be between -180 and 180"

        # Use specified resolution or VGRID instance resolution
        target_resolution = resolution if resolution is not None else self.resolution

        # Calculate cell size for the target resolution
        # Base cell size is the VGRID instance's cell size
        base_cell_size = self.cell_size
        if target_resolution >= 0:
            # For positive resolutions, divide both longitude and latitude cell sizes by sqrt(aperture)^resolution
            # This ensures total cells increase by aperture factor, not aperture^2
            factor = (self.aperture**0.5) ** target_resolution
            target_cell_size_lon = base_cell_size / factor
            target_cell_size_lat = base_cell_size / factor
        else:
            # For negative resolutions, multiply both longitude and latitude cell sizes by sqrt(aperture)^|resolution|
            factor = (self.aperture**0.5) ** abs(target_resolution)
            target_cell_size_lon = base_cell_size * factor
            target_cell_size_lat = base_cell_size * factor
        target_total_columns = int(360 / target_cell_size_lon)
        target_total_rows = int(180 / target_cell_size_lat)

        # If crossing the anti-meridian, split it up and combine
        if left > right:
            east = self.tiles_for_bounding_box(
                left, bottom, 180.0, top, target_resolution
            )
            west = self.tiles_for_bounding_box(
                -180.0, bottom, right, top, target_resolution
            )
            return east + west

        # Move coordinates so we can compute percentages
        left += 180
        right += 180
        bottom += 90
        top += 90

        tiles = []
        # For each column (longitude)
        for x in range(
            int(left / target_cell_size_lon), int(right / target_cell_size_lon) + 1
        ):
            # For each row (latitude)
            for y in range(
                int(bottom / target_cell_size_lat), int(top / target_cell_size_lat) + 1
            ):
                # Ensure we don't exceed grid bounds
                if (
                    x >= 0
                    and x < target_total_columns
                    and y >= 0
                    and y < target_total_rows
                ):
                    # Calculate tile index based on the target resolution's grid
                    tile_index = y * target_total_columns + x
                    tiles.append((target_resolution, tile_index))

        return tiles

    def __repr__(self) -> str:
        return f"VGRID(cell_size={self.cell_size}, resolution={self.resolution}, aperture={self.aperture})"

    def __eq__(self, other) -> bool:
        if not isinstance(other, VGRID):
            return False
        return (
            self.cell_size == other.cell_size
            and self.resolution == other.resolution
            and self.aperture == other.aperture
        )
