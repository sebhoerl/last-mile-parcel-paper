import attrs
import xarray as xr
import numpy as np

@attrs.define
class VehicleType:
    distance_cost_EUR_per_m: float = attrs.field(
        validator = attrs.validators.ge(0.0))

    vehicle_cost_EUR_per_day: float = attrs.field(
        validator = attrs.validators.ge(0.0))
    
    capacity: int = attrs.field(
        validator = attrs.validators.gt(0))
    
    maximum_speed_km_h: float = attrs.field(
        validator = attrs.validators.gt(0.0))

    active_time_h: float = attrs.field(
        validator = attrs.validators.ge(0.0))
    
    network: str = "default"

    maximum_travel_time_h: float = attrs.field(
        validator = attrs.validators.ge(0.0), default = np.inf)
    
    maximum_distance_km: float = attrs.field(
        validator = attrs.validators.ge(0.0), default = np.inf)

    maximum_travel_time_per_tour_h: float = attrs.field(
        validator = attrs.validators.ge(0.0), default = np.inf)
    
    maximum_distance_per_tour_km: float = attrs.field(
        validator = attrs.validators.ge(0.0), default = np.inf)
    
    maximum_vehicles: int = None

@attrs.define
class Shipment:
    id: str
    node: int

@attrs.define
class Skills:
    available_skills: dict[int, list[str]] = attrs.field()
    node_skills: list[list[int]] = attrs.field()
    network_skills: dict[str, list[int]] = attrs.field()
    shipment_skills: list[list[int]] = attrs.field()

@attrs.define
class Problem:
    vehicle_types: dict[str, VehicleType] = attrs.field(
        validator = attrs.validators.min_len(1))

    matrices: dict[str, xr.DataArray] = attrs.field(
        validator = attrs.validators.min_len(1))

    pickup_duration: float = attrs.field(
        validator = attrs.validators.ge(0.0))
    
    delivery_duration: float = attrs.field(
        validator = attrs.validators.ge(0.0))
    
    shipments: list[Shipment] = attrs.field(
        validator = attrs.validators.min_len(1))
    
    depot_node: int = attrs.field(
        validator = attrs.validators.ge(0.0))
    
    travel_time_factor: float = attrs.field(default = 1.0)

    use_tour_constraints: bool = attrs.field(default = True)

    _skills: Skills = attrs.field(kw_only = True, default = None)

    @matrices.validator
    def validate_matrices(self, attribute, value):
        required_matrices = set([vt.network for vt in self.vehicle_types.values()])

        if len(required_matrices - self.matrices.keys()) > 0:
            raise ValueError("Missing matrices: {}".format(required_matrices - self.matrices.keys()))
        
        for name, matrix in self.matrices.items():
            if matrix.dims != ("origin", "destination", "attribute"):
                raise ValueError("Wrong matrix format for {}".format(name))
            
            if set(matrix.coords["attribute"].values) != set(["distance", "travel_time"]):
                raise ValueError("Matrix does not containt distance and travel time: {}".format(name))

            node_sequence = matrix.coords["origin"].values
            
            if not np.all(node_sequence == matrix.coords["origin"].values):
                raise ValueError("Not a proper square matrix for {}".format(name))
            
            if not np.all(node_sequence == matrix.coords["destination"].values):
                raise ValueError("Not a proper square matrix for {}".format(name))

    @shipments.validator
    def validate_shipments(self, attribute, value):
        nodes = set([shipment.node for shipment in self.shipments])

        for matrix in self.matrices.values():
            nodes -= set(matrix.coords["origin"].values)

        if len(nodes) > 0:
            raise ValueError("Some nodes are not covered by any matrix: {}".format(nodes))
            
    @depot_node.validator
    def validate_depot_node(self, attributes, value):
        for matrix in self.matrices.values():
            if self.depot_node in matrix.coords["origin"]:
                return
        
        raise ValueError("Depot node not covered by any matrix")

    def _reindex_matrix(self, matrix):
        node_sequence = self.get_node_sequence()
        fill_value = self._matrix_forbidden_value(matrix)

        #if len(matrix.coords["origin"]) > 0:
        return matrix.reindex({
            "origin": node_sequence,
            "destination": node_sequence
        }, fill_value = fill_value)

    def _matrix_forbidden_value(self, matrix):
        return matrix.values.max() * 10.0 if len(matrix.coords["origin"]) > 0 else 1000.0

    def get_distance_matrix(self, vehicle_type):
        network = self.vehicle_types[vehicle_type].network
        return self._reindex_matrix(self.matrices[network].sel(attribute = "distance", drop = True))
    
    def get_cost_matrix(self, vehicle_type):
        cost_per_distance = self.vehicle_types[vehicle_type].distance_cost_EUR_per_m
        return self.get_distance_matrix(vehicle_type) * cost_per_distance
    
    def get_travel_time_matrix(self, vehicle_type):
        vehicle_type = self.vehicle_types[vehicle_type]
        matrix = self.matrices[vehicle_type.network]

        distances = matrix.sel(attribute = "distance", drop = True)
        vehicle_travel_times = 3.6 * distances / vehicle_type.maximum_speed_km_h

        travel_times = matrix.sel(attribute = "travel_time", drop = True)
        travel_times = np.maximum(travel_times, vehicle_travel_times)
        travel_times *= self.travel_time_factor

        return self._reindex_matrix(travel_times)
    
    def get_vehicle_unit_cost(self, vehicle_type):
        vehicle_type = self.vehicle_types[vehicle_type]
        return vehicle_type.vehicle_cost_EUR_per_day
    
    def get_node_sequence(self):
        all_nodes = set()

        for matrix in self.matrices.values():
            all_nodes |= set(matrix.coords["origin"].values)

        return list(sorted(list(all_nodes))) 
    
    def get_network_sequence(self):
        return list(sorted(self.matrices.keys()))
    
    def get_skills(self):
        if self._skills is not None:
            return self._skills

        # Obtain all nodes
        node_sequence = self.get_node_sequence()
        networks = self.get_network_sequence()

        # Matrix nodes
        network_nodes = {
            network: set(self.matrices[network].coords["origin"].values)
            for network in networks
        }

        # Obtain skills as combinations of networks
        all_skills = []
        node_skills = {}

        for node in node_sequence:
            list_skill = []

            for network in networks:
                if node in network_nodes[network]:
                    list_skill.append(network)

            if len(list_skill) < len(networks):
                list_skill = tuple(list_skill)

                if not list_skill in all_skills:
                    all_skills.append(list_skill)

                node_skills[node] = [all_skills.index(list_skill)]
            else:
                node_skills[node] = []

        # Obtain vehicle skills depending on membership in network skills
        network_skills = {}

        for network in networks:
            current_skills = []

            for skill, list_skill in enumerate(all_skills):
                if network in list_skill:
                    current_skills.append(skill)
            
            network_skills[network] = current_skills

        available_skills = {
            skill: list_skill
            for skill, list_skill in enumerate(all_skills)
        }

        # Obtain shipment skills depending on the depot and destination node
        depot_skills = set(node_skills[self.depot_node])
        shipment_skills = []

        for shipment in self.shipments:
            local_skills = set(node_skills[shipment.node]) | depot_skills
            shipment_skills.append(list(sorted(list(local_skills))))

        self._skills = Skills(available_skills, node_skills, network_skills, shipment_skills)
        return self._skills

    def verify_sequence(self, network, nodes):
        return len(set(nodes) - set(self.matrices[network].coords["origin"].values)) == 0
