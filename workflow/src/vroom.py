import solver
import numpy as np
import json, orjson
import logging
import packaging.version
import hashlib
import tempfile
import subprocess as sp
import os
import attrs
import time
import pandas as pd

logger = logging.getLogger(__name__)

@attrs.define
class VroomInstance:
    ORJSON_OPTIONS = orjson.OPT_SORT_KEYS | orjson.OPT_SERIALIZE_NUMPY

    problem: solver.Problem = attrs.field()
    seed: int = attrs.field(default = 0)
    cost_factor: float = attrs.field(default = 1000.0)
    _data: map = attrs.field()

    @_data.default
    def _initialize_data(self):
        data = {}

        # Matrices
        matrices = {}
        data["matrices"] = matrices

        need_distance_matrices = False
        for vehicle_type in self.problem.vehicle_types.values():
            if np.isfinite(vehicle_type.maximum_distance_per_tour_km):
                need_distance_matrices = True

            if np.isfinite(vehicle_type.maximum_distance_km):
                need_distance_matrices = True

        for vehicle_type in sorted(list(self.problem.vehicle_types.keys())):
            cost_matrix = self.problem.get_cost_matrix(vehicle_type)
            travel_time_matrix = self.problem.get_travel_time_matrix(vehicle_type)

            matrices[vehicle_type] = {
                "durations": (travel_time_matrix.values).astype(int),
                "costs": (cost_matrix.values * self.cost_factor).astype(int),
            }

            if need_distance_matrices:
                matrices[vehicle_type]["distances"] = self.problem.get_distance_matrix(vehicle_type).values.astype(int)

        # Skills
        skills = self.problem.get_skills()

        # Shipments
        shipments = []
        data["shipments"] = shipments

        node_sequence = self.problem.get_node_sequence()
        # is_depot_electric = self.problem.depot_node in self.problem.electric_nodes

        for shipment_index, shipment in enumerate(self.problem.shipments):
            # is_electric = is_depot_electric or (shipment.node in self.problem.electric_nodes)

            shipments.append({
                "pickup": {
                    "id": shipment_index,
                    "location_index": node_sequence.index(self.problem.depot_node),
                    "service": int(self.problem.pickup_duration)
                },
                "delivery": {
                    "id": shipment_index,
                    "location_index": node_sequence.index(shipment.node),
                    "service": int(self.problem.delivery_duration),
                },
                "amount": [1],
                "skills": skills.shipment_skills[shipment_index]
            })

        # Vehicles
        data["vehicles"] = []

        return data

    def set_vehicles(self, vehicles_per_type: dict[str, int]):
        vehicles = []

        node_sequence = self.problem.get_node_sequence()
        depot_index = node_sequence.index(self.problem.depot_node)
        skills = self.problem.get_skills()

        for vehicle_type, count in vehicles_per_type.items():
            daily_cost = self.problem.get_vehicle_unit_cost(vehicle_type)
            options = self.problem.vehicle_types[vehicle_type]

            for k in range(count):
                vehicle = {
                    "profile": vehicle_type,
                    "start_index": depot_index,
                    "end_index": depot_index,
                    "capacity": [options.capacity],
                    "time_window": [0, int(options.active_time_h * 3600)],
                    "costs": {
                        "fixed": int(self.cost_factor * daily_cost)
                    },
                    "skills": skills.network_skills[options.network]
                }

                if np.isfinite(options.maximum_travel_time_h):
                    vehicle["max_travel_time"] = int(options.maximum_travel_time_h * 3600)

                if np.isfinite(options.maximum_distance_km):
                    vehicle["max_distance"] = int(options.maximum_distance_km * 1000)

                if self.problem.use_tour_constraints:
                    # Custom addition to VROOM
                    # https://github.com/VROOM-Project/vroom/issues/918
                
                    if np.isfinite(options.maximum_travel_time_per_tour_h):
                        vehicle["max_tour_travel_time"] = int(options.maximum_travel_time_per_tour_h * 3600)
                    
                    if np.isfinite(options.maximum_distance_per_tour_km):
                        vehicle["max_tour_distance"] = int(options.maximum_distance_per_tour_km * 1000)

                vehicles.append(vehicle)

        # Randomize IDs
        random = np.random.RandomState(self.seed)
        random.shuffle(vehicles)

        for k in range(len(vehicles)):
            vehicles[k]["id"] = k

        self._data["vehicles"] = vehicles

    def write(self, path):
        with open(path, "wb+") as f:
            f.write(orjson.dumps(self._data, option = self.ORJSON_OPTIONS))

    def hash(self):
        return hashlib.md5(orjson.dumps(self._data, option = self.ORJSON_OPTIONS)).hexdigest()

    def create_solution(self, solution: dict):
        vehicles_by_id = {
            vehicle["id"]: vehicle
            for vehicle in self._data["vehicles"]
        }

        df_stops = []

        for route in solution["routes"]:
            vehicle_type = vehicles_by_id[route["vehicle"]]["profile"]

            location_indices = []
            stop_types = []

            arrival_times = []
            departure_times = []

            shipment_indices = []

            current_time = 0.0
            current_travel_time = 0.0

            vehicle_cost = vehicles_by_id[route["vehicle"]]["costs"]["fixed"]
            current_cost = vehicle_cost

            for step in route["steps"]:
                # track movement
                if len(location_indices) > 0:
                    previous_location_index = location_indices[-1]
                    current_location_index = step["location_index"]

                    step_travel_time = self._data["matrices"][vehicle_type]["durations"][
                        previous_location_index, current_location_index]
                    
                    current_travel_time += step_travel_time
                    current_time += step_travel_time

                    current_cost += self._data["matrices"][vehicle_type]["costs"][
                        previous_location_index, current_location_index]
                    
                # track location
                location_indices.append(step["location_index"])
                arrival_times.append(current_time)

                # verify timing
                assert step["duration"] == current_travel_time
                assert step["arrival"] == current_time

                # service timing
                current_time += step["service"] if "service" in step else 0.0
                departure_times.append(current_time)

                # additional information
                stop_types.append(step["type"])
                shipment_indices.append(step["id"] if "id" in step else -1)
            
            # verify cost
            assert current_cost == route["cost"]

            df_partial = pd.DataFrame({
                "stop_sequence": list(range(len(location_indices))),
                "location_index": location_indices,
                "stop_type": stop_types,
                "arrival_time": arrival_times,
                "departure_time": departure_times,
                "shipment_index": shipment_indices
            })

            df_partial["vehicle_id"] = route["vehicle"]
            df_partial["vehicle_type"] = vehicle_type

            df_stops.append(df_partial)

        instance_information = {
            "feasible": solution["feasible"],
            "computing_times": solution["summary"]["computing_times"],
            "retries": solution["retries"]
        }

        return pd.concat(df_stops), instance_information

@attrs.define
class VroomExecutable:
    REQUIRED_VERSION = "1.12.0"

    executable: str = attrs.field(default = "vroom")
    working_directory: str = attrs.field(default = None)
    time_limit_s: int = attrs.field(default = None)
    threads: int = attrs.field(default = None)
    exploration: int = attrs.field(default = 5)

    def run(self, instance: VroomInstance):
        # Verfiy VROOM
        valid_version = False

        version_output = sp.check_output([
            self.executable, "-h"
        ]).decode("utf-8")

        for line in version_output.split("\n"):
            if line .startswith("Version: "):
                version = line.replace("Version: ", "")
                valid_version = packaging.version.parse(version) >= packaging.version.parse(self.REQUIRED_VERSION)

        if not valid_version:
            raise RuntimeError("Wrong version of VROOM")
        
        # Prepare directory
        hash = instance.hash()

        working_directory = self.working_directory
        temporary_directory = None 

        if self.working_directory is None:
            temporary_directory = tempfile.TemporaryDirectory()
            working_directory = temporary_directory.name

        problem_path = "{}/vroom_{}_problem.json".format(working_directory, hash)
        solution_path = "{}/vroom_{}_solution.json".format(working_directory, hash)

        # Write problem
        writing_start_time = time.time()
        instance.write(problem_path)
        writing_end_time = time.time()

        # Prepare command
        command = [
            self.executable,
            "-i", problem_path,
            "-o", solution_path,
            "-x", str(self.exploration)
        ]

        if not self.time_limit_s is None:
            command += ["-l", str(self.time_limit_s)]

        if not self.threads is None:
            command += ["-t", str(self.threads)]

        # Execute command
        try:
            sp.check_call(command)
        except sp.CalledProcessError:
            os.remove(problem_path)

            if temporary_directory is not None:
                temporary_directory.cleanup()
            
            return { "feasible": False }

        reading_start_time = time.time()
        with open(solution_path, "r") as f:
            solution = json.load(f)
        reading_end_time = time.time()

        os.remove(problem_path)
        os.remove(solution_path)
        
        if temporary_directory is not None:
            temporary_directory.cleanup()

        solution["summary"]["computing_times"]["python_writing"] = (writing_end_time - writing_start_time) * 1e3
        solution["summary"]["computing_times"]["python_reading"] = (reading_end_time - reading_start_time) * 1e3
        solution["feasible"] = True
        return solution

@attrs.define
class VroomSolver:
    instance: VroomInstance = attrs.field()
    executable: VroomExecutable = attrs.field()
    seed: int = attrs.field(default = 0)
    vehicles_per_type: int = attrs.field(default = 40)
    maximum_retries: int = attrs.field(default = 5)

    def solve(self):
        maximum_vehicles_per_type = {
            vehicle_type: np.inf if options.maximum_vehicles is None else options.maximum_vehicles
            for vehicle_type, options in self.instance.problem.vehicle_types.items()
        }

        vehicles_per_type = {
            vehicle_type: min(self.vehicles_per_type, maximum_vehicles_per_type[vehicle_type])
            for vehicle_type in self.instance.problem.vehicle_types.keys()
        }

        retries = 0
        solution = None
        continue_solving = True

        while continue_solving:
            self.instance.set_vehicles(vehicles_per_type)

            print("Solving instance with {}".format(vehicles_per_type))
            solution = self.executable.run(self.instance)

            increase_vehicle_types = set()

            # solution not feasible, try to increase all vehicle counts
            if not self._is_feasible(solution):
                print("No feasible solution found for instance")
                
                for vehicle_type in vehicles_per_type.keys():
                    increase_vehicle_types.add(vehicle_type)
          
            # solution feasible, check if any vehicle type reached the count
            else:
                used_vehicles_per_type = self._get_vehicles_per_type(self.instance, solution)
                print("Feasible solution for isntance: {}".format(used_vehicles_per_type))
                
                for vehicle_type in vehicles_per_type.keys():
                    if used_vehicles_per_type[vehicle_type] == vehicles_per_type[vehicle_type]:
                        # for this vehicle_type we reached the indicated bound, try to increase
                        increase_vehicle_types.add(vehicle_type)

            if len(increase_vehicle_types) == 0:
                # found a solution and no need to increase vehicles
                continue_solving = False

            else:
                continue_solving = False # by default, abort

                for vehicle_type in increase_vehicle_types:
                    if vehicles_per_type[vehicle_type] < maximum_vehicles_per_type[vehicle_type]:
                        # we could increase some vehicle count
                        vehicles_per_type[vehicle_type] = min(
                            vehicles_per_type[vehicle_type] + self.vehicles_per_type, maximum_vehicles_per_type[vehicle_type])
                        continue_solving = True # continue

            retries += 1

            if retries > self.maximum_retries:
                continue_solving = False

        if not self._is_feasible(solution):
            solution = None
        else:
            solution["retries"] = retries
            solution = self.instance.create_solution(solution)

        return solution

    def _is_feasible(self, solution):
        return len(solution["unassigned"]) == 0 and solution["feasible"]

    def _get_vehicles_per_type(self, instance, solution):
        vehicle_types = {
            vehicle["id"]: vehicle["profile"]
            for vehicle in instance._data["vehicles"]
        }

        counts = { vehicle_type: 0 for vehicle_type in instance._data["matrices"].keys() } 
        for route in solution["routes"]:
            vehicle_type = vehicle_types[route["vehicle"]]
            counts[vehicle_type] = counts[vehicle_type] + 1

        return counts
