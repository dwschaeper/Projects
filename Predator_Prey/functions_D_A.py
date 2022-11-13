import random as r
import matplotlib.pyplot as plt
from math import dist, sqrt


# ----------------------------------------------------------------------------------------------------------------------
def initialize_birds(canvas, attributes, hawk_quant, dove_quant, starting_bird_energy, width, height):
    """=================================================================================================================
    This function creates the initial birds to be competing for food, and creates a dictionary with all of their
    properties.

    :param canvas: the canvas the birds are placed in
    :param attributes: empty dictionary to have attributes placed in
    :param hawk_quant: number of simulation hawks
    :param dove_quant: number of simulation doves
    :param starting_bird_energy: the amount of energy a bird starts withS
    :param width: the width of the canvas
    :param height: the height of the canvas
    :return: dictionary, the bird attributes
    ================================================================================================================="""
    for i in range(0, dove_quant):
        # create a numbered tag for each individual bird
        dove_id = 'dove' + str(i)

        # randomly generate the top left corner coordinate
        dove_x = r.randint(100, width - 100)
        dove_y = r.randint(100, height - 100)

        # create the bird object
        canvas.create_oval(dove_x, dove_y, dove_x + 30, dove_y + 30, fill='blue', outline='black',
                           tag=('bird', dove_id))

        # NEW
        attributes[dove_id] = {'alive': 'Yes', 'energy': starting_bird_energy, 'location': [dove_x + 15, dove_y - 15],
                               'vector': None,
                               'type': 'Dove', 'reproduced': 150, 'time_eating': 0, 'eating?': False,
                               'food_coordinate?': None,
                               'food_name': None, 'food_target': None}

    for i in range(0, hawk_quant):
        # create a numbered tag for each individual bird
        hawk_id = 'hawk' + str(i)

        # randomly generate the top left corner coordinate
        hawk_x = r.randint(100, width - 100)
        hawk_y = r.randint(100, height - 100)

        # create the bird object
        canvas.create_oval(hawk_x, hawk_y, hawk_x + 30, hawk_y + 30, fill='red', outline='black', tag=('bird', hawk_id))

        # NEW
        attributes[hawk_id] = {'alive': 'Yes', 'energy': starting_bird_energy, 'location': [hawk_x + 15, hawk_y - 15],
                               'vector': None,
                               'type': 'Hawk', 'reproduced': 150, 'time_eating': 0, 'eating?': False,
                               'food_coordinate?': None,
                               'food_name': None, 'food_target': None}

    return attributes


# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
def initialize_food(canvas, attributes, food_number, width, height):
    """=================================================================================================================
    This function creates the initial food that the birds are competing for, and creates a dictionary of the properties
    of the food objects.

    :param canvas: the canvas the food is placed in
    :param attributes: an empty dictionary to place the attributes in
    :param food_number: the amount of food items to be in the sim
    :param width: the width of the canvas
    :param height: the height of the canvas
    :return: dictionary, the food attributes
    ================================================================================================================="""
    for i in range(0, food_number):
        # create a numbered id for individual food pieces
        food_id = 'food' + str(i)

        # randomly generate the top left corner of the food object
        x = r.randint(100, width - 100)
        y = r.randint(100, height - 100)

        # create the food object
        canvas.create_oval(x, y, x + 10, y + 10, fill='green', outline='black', tag=('food', food_id))

        attributes[food_id] = {'active': 'yes', 'location': [x, y], 'being_eaten': 'No'}
    return attributes


# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
def collect_coords(canvas):
    """=================================================================================================================
    UNUSED
    This function collects the coords of the food and bird objects.

    :param canvas: the canvas object the food and birds are placed in
    :return: food_coords, hawk_coords, dove_coords: lists of the coords
    ================================================================================================================="""
    # initialize variables
    i = 0
    food_coords = []
    hawk_coords = {}
    dove_coords = {}

    # loop to get all of the coordinates for the food objects. Multiple loops because through reproduction and death,
    # the total number of objects for each will be variable
    while True:
        tag = 'food' + str(i)
        coord = canvas.coords(tag)[0:2]
        if coord:
            food_coords.append([tag, coord])
            i += 1
        else:
            i = 0
            break

    # loop to get the coordinates for the hawk objects
    while True:
        tag = 'hawk' + str(i)
        coord = canvas.coords(tag)[0:2]
        if coord:
            # hawk_coords.append([tag,coord])
            hawk_coords[tag] = coord
            i += 1
        else:
            i = 0
            break

    # loop to get the coordinates for the dove object
    while True:
        tag = 'dove' + str(i)
        coord = canvas.coords(tag)[0:2]
        if coord:
            # dove_coords.append([tag,coord])
            dove_coords[tag] = coord
            i += 1
        else:
            break

    return food_coords, hawk_coords, dove_coords


# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
def calculate_vector(bird_attributes, food_attributes):
    """=================================================================================================================
    Calculates the distance between the bird and the food to find the food object that is closest the bird. The food
    object that is closest to the bird will be the target for the direction vector that moves the bird toward the food.

    :param bird_attributes: a dictionary of all the attributes of the current birds
    :param food_attributes:a dictionary of all the attributes of the current food
    :return: An updated bird_attributes with unit vectors
    ================================================================================================================="""
    for bird in bird_attributes:
        min_distance = 6000
        # get the coords of the birds
        bird_coords = bird_attributes[bird]['location']
        for food in food_attributes:
            # get coords of the food

            # take_food=bool(r.getrandbits(1))

            if food_attributes[food]['active'] == 'No':
                continue

            # elif food_attributes[food]['being_eaten'] == 'Yes':
            #    continue

            else:
                food_coords = food_attributes[food]['location']
                distance = dist(bird_coords, food_coords)
                # find the closest food to the bird
                if distance < min_distance:
                    min_distance = distance
                    vector = [food_coords[0] - bird_coords[0], food_coords[1] - bird_coords[1]]
                    magnitude = sqrt((vector[0] ** 2) + (vector[1] ** 2))
                    unit_vector_1 = [vector[0] / magnitude, vector[1] / magnitude]
                    # set the movement vector the unit vector between the bird and the closest food
                    bird_attributes[bird]['vector'] = unit_vector_1
    return bird_attributes


# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
def reproduction(canvas, x, y, bird_attributes, species, energy, mother_tag):
    """=================================================================================================================
    This function is triggered when the bird is well-fed as defined by the user.

    :param canvas: The canvas object the sim is in
    :param x: x coord of parent bird
    :param y: y coord of parent bird
    :param bird_attributes: list of birds and their attributes
    :param species: the type of bird being reproduced
    :param energy: the energy the bird starts with
    :param mother_tag: the parent of the progeny

    :return: NULL
    ================================================================================================================="""
    progeny_x = r.randint(int(x) - 50, int(x) + 50)
    progeny_y = r.randint(int(y) - 50, int(y) + 50)

    if species == 'Hawk':
        tag = 'hawk' + str(len(bird_attributes))
        canvas.create_oval(progeny_x, progeny_y, progeny_x + 30, progeny_y + 30, fill='red', outline='black',
                           tag=('bird', tag))
    else:
        tag = 'dove' + str(len(bird_attributes))
        canvas.create_oval(progeny_x, progeny_y, progeny_x + 30, progeny_y + 30, fill='blue', outline='black',
                           tag=('bird', tag))
    bird_attributes[tag] = {'alive': 'Yes', 'energy': energy, 'location': [progeny_x + 15, progeny_y + 15],
                            'vector': [-bird_attributes[mother_tag]['vector'][0],
                                       -bird_attributes[mother_tag]['vector'][1]],
                            'type': species, 'reproduced': 150, 'time_eating': 0, 'eating?': False,
                            'food_coordinate?': None,
                            'food_name': None, 'food_target': None}


# ----------------------------------------------------------------------------------------------------------------------


def simulation_visualization(bird_attributes_list):
    # hawk_dove_ratios = [[ii[:-1:] for ii in i] for i in bird_attributes_list]
    simulation_cycles = 0
    hawk_dove_ratio_cycle = []
    for i in bird_attributes_list:
        simulation_cycles += 1
        hawk_count = 0
        dove_count = 0
        for ii in i:
            if ii == 'hawk':
                hawk_count += 1
            elif ii == 'dove':
                dove_count += 1
        hawk_dove_ratio_cycle.append(hawk_count / (hawk_count + dove_count))

    total_cycles = list(range(0, simulation_cycles))

    # print(hawk_dove_ratio_cycle)
    plt.plot(total_cycles, hawk_dove_ratio_cycle, 'r')
    plt.xlabel('number of cycles')
    plt.ylabel('hawk_count / (hawk_count + dove_count), (0-1)')
    plt.ylim(0, 1)
    plt.show()
