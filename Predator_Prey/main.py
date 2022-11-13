import tkinter as tk
import random as r
import functions_D_A


# ----------------------------------------------------------------------------------------------------------------------


def run_sim(canvas, root, bird_attributes, food_attributes, energy_decrease, time_step, food_calories,
            reproduction_energy, width, height, eating_time, bird_attribute_list):
    """=================================================================================================================
    Runs the simulation

    :param canvas: the canvas object
    :param root: the main window
    :param bird_attributes: the dictionary of birds
    :param food_attributes: the dictionary of food
    :param energy_decrease: the decrement amount of energy for moving
    :param time_step: the time between each tic
    :param food_calories: the amount of energy gained for eating
    :param reproduction_energy: the amount of energy required to reproduce
    :param width: the width of the canvas object
    :param height:the height of the canvas object
    :param eating_time: the time it takes for a bird to eat

    :return: null
    ================================================================================================================="""

    functions_D_A.calculate_vector(bird_attributes, food_attributes)
    birds = []
    # bird_attribute_list.append(bird_attributes)

    for j in list(bird_attributes):
        # if the bird is dead, remove it

        if bird_attributes[j]['alive'] == 'No':
            canvas.delete(j)
        # continue with the sim if the bird is not dead
        else:
            vector = bird_attributes[j]['vector']

            # move the bird by its unit vector if the bird is not currently eating

            if bird_attributes[j]['time_eating'] == 0 and not bird_attributes[j]['eating?']:
                canvas.move(j, vector[0], vector[1])

            # update the bird coordinates in the attributes
            bird_coords = canvas.coords(j)[0:2]
            bird_attributes[j]['location'] = bird_coords

            # iterate through the food items to see if a bird should start eating
            for ii in list(food_attributes):
                # need to check if the food item is active or an eaten on still in the list, explains why its still
                # there later
                if food_attributes[ii]['active'] == 'No':
                    continue
                else:
                    # check if the bird is already eating, don't want to reset it each time
                    if bird_attributes[j]['eating?']:
                        continue
                    # if the food object is being eaten, don't let other birds eat it as well.
                    # this may need to change depending on what you plan for interactions

                    # elif food_attributes[ii]['being_eaten'] == 'Yes':
                    #    continue

                    # better way to write to check if an object is touching with just one if statement, i knew there was
                    # a way, just took a bit to remember
                    elif int(food_attributes[ii]['location'][0]) - 15 <= canvas.coords(j)[0:2][0] <= int(
                            food_attributes[ii]['location'][0]) + 15 and int(food_attributes[ii]['location'][1]) - 15 \
                            <= canvas.coords(j)[0:2][1] <= int(food_attributes[ii]['location'][1]) + 15 \
                            and not food_attributes[ii]['being_eaten'] == 'Yes':
                        # get the location of the food the bird is eating
                        bird_attributes[j]['food_coordinate?'] = [food_attributes[ii]['location'][0],
                                                                  food_attributes[ii]['location'][1]]
                        bird_attributes[j]['eating?'] = True
                        bird_attributes[j]['food_name'] = ii
                        food_attributes[ii]['being_eaten'] = 'Yes'

                    # Signals that the bird has found food that another bird was devouring
                    elif int(food_attributes[ii]['location'][0]) - 15 <= canvas.coords(j)[0:2][0] <= int(
                            food_attributes[ii]['location'][0]) + 15 and int(food_attributes[ii]['location'][1]) - 15 \
                            <= canvas.coords(j)[0:2][1] <= int(food_attributes[ii]['location'][1]) + 15 \
                            and food_attributes[ii]['being_eaten'] == 'Yes':

                        # Gets the location of the food the bird is eating
                        bird_attributes[j]['food_coordinate?'] = [food_attributes[ii]['location'][0],
                                                                  food_attributes[ii]['location'][1]]

                        # Work in dove-dove interaction

                        for bird in list(bird_attributes):
                            if bird != j:
                                if bird_attributes[j]['food_coordinate?'] \
                                        == bird_attributes[bird]['food_coordinate?']:

                                    bird_attributes[j]['food_name'] = ii
                                    bird_attributes[bird]['food_name'] = ii

                                    bird_attributes[j]['eating?'] = False
                                    bird_attributes[j]['time_eating'] = 0

                                    bird_attributes[bird]['eating?'] = False
                                    bird_attributes[bird]['time_eating'] = 0

                                    if bird_attributes[j]['type'] == 'Dove' and bird_attributes[bird]['type'] == 'Dove':
                                        bird_attributes[j]['energy'] += food_calories / 2
                                        bird_attributes[bird]['energy'] += food_calories / 2
                                    elif bird_attributes[j]['type'] == 'Hawk' and bird_attributes[bird][
                                        'type'] == 'Dove':
                                        bird_attributes[j]['energy'] += food_calories
                                    elif bird_attributes[j]['type'] == 'Dove' and bird_attributes[bird][
                                        'type'] == 'Hawk':
                                        bird_attributes[bird]['energy'] += food_calories
                                    elif bird_attributes[j]['type'] == 'Hawk' and bird_attributes[bird][
                                        'type'] == 'Hawk':
                                        pass
                                    # get rid of the eaten food

                                    canvas.delete(bird_attributes[j]['food_name'])
                                    food_attributes[bird_attributes[j]['food_name']]['active'] = 'No'

                                    canvas.delete(bird_attributes[bird]['food_name'])
                                    food_attributes[bird_attributes[bird]['food_name']]['active'] = 'No'

                                    bird_attributes[j]['food_coordinate?'] = None
                                    bird_attributes[bird]['food_coordinate?'] = None
                                    # generate the replacement food
                                    x = r.randint(50, width - 50)
                                    y = r.randint(50, height - 50)
                                    # need a new tag because deletion won't properly with 2 items with the same tag

                                    # could be problematic

                                    tag = 'food' + str(len(food_attributes))
                                    canvas.create_oval(x, y, x + 10, y + 10, fill='green', outline='black',
                                                       tag=('food', tag))
                                    food_attributes[tag] = {'active': 'yes', 'location': [x, y], 'being_eaten': 'No'}
                                    break

            # check if a bird can reproduce
            if bird_attributes[j]['energy'] > reproduction_energy:
                bird_attributes[j]['energy'] = (bird_attributes[j]['energy']) / 2
                bird_attributes[j]['reproduced'] = 150
                # generate the new bird
                functions_D_A.reproduction(canvas=canvas, x=bird_coords[0], y=bird_coords[1],
                                           bird_attributes=bird_attributes,
                                           species=bird_attributes[j]['type'],
                                           energy=(bird_attributes[j]['energy']) / 2,
                                           mother_tag=j)

            # decrease birds energy based on movement cost
            bird_attributes[j]['energy'] = bird_attributes[j]['energy'] - energy_decrease
            # kill bird if it runs out of energy
            if bird_attributes[j]['energy'] <= 0:
                bird_attributes[j]['alive'] = 'No'

            # if the bird is eating, increment the eating duration
            # when the duration is reached, create a new food item, give food energy, and destroy the original food item
            if bird_attributes[j]['eating?']:
                bird_attributes[j]['time_eating'] += 1
                # if the food is all the way eaten
                if bird_attributes[j]['time_eating'] == eating_time:
                    bird_attributes[j]['eating?'] = False
                    bird_attributes[j]['time_eating'] = 0
                    bird_attributes[j]['energy'] += food_calories
                    # get rid of the eaten food
                    canvas.delete(bird_attributes[j]['food_name'])
                    food_attributes[bird_attributes[j]['food_name']]['active'] = 'No'
                    # generate the replacement food
                    x = r.randint(50, width - 50)
                    y = r.randint(50, height - 50)
                    # need a new tag because deletion won't properly with 2 items with the same tag
                    tag = 'food' + str(len(food_attributes))
                    canvas.create_oval(x, y, x + 10, y + 10, fill='green', outline='black', tag=('food', tag))
                    food_attributes[tag] = {'active': 'yes', 'location': [x, y], 'being_eaten': 'No'}
        birds.append(j[0:4])
    bird_attribute_list.append(birds)
    # print(len(bird_attribute_list))

    root.after(time_step, run_sim, canvas, root, bird_attributes, food_attributes, energy_decrease, time_step,
               food_calories, reproduction_energy, width, height, eating_time, bird_attribute_list)


# ----------------------------------------------------------------------------------------------------------------------


# ----------------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    """=================================================================================================================
    This is the main program.

    :return The GUI object
    ================================================================================================================="""

    # global bird_attributes_list
    bird_attributes_list = []

    # Simulation parameters
    # field dimensions
    width = 1400
    height = 1000
    # Number of doves and hawks respectively
    dove_quant = 3
    hawk_quant = 20
    # Number of food particles to remain in the simulation at all times
    food_number = 25
    # Starting energy for each hawk and dove
    starting_bird_energy = 400
    # Reproduction energy
    reproduction_energy = 600
    # decrease in energy per step
    energy_decrease = 0.1
    # time between steps, simulation will slow when many creatures are inside regardless
    time_step = 10
    # Determines the number of calories per food
    food_calories = 100
    # Amount of time spent eating
    eating_time = 500

    root = tk.Tk()
    root.geometry(str(width) + 'x' + str(height))
    root.title('Hawks and Doves Sim')
    # create canvas
    canvas = tk.Canvas(root, width=width, height=height)
    canvas.pack()
    # initialize holder for attributes
    bird_attributes = {}
    food_attributes = {}
    bird_attributes = functions_D_A.initialize_birds(canvas, bird_attributes, hawk_quant, dove_quant,
                                                     starting_bird_energy, width, height)
    food_attributes = functions_D_A.initialize_food(canvas, food_attributes, food_number, width, height)
    # food_coords, hawk_coords, dove_coords = Functions_D.collect_coords(canvas)
    bird_attributes = functions_D_A.calculate_vector(bird_attributes, food_attributes)

    run_sim(canvas, root, bird_attributes, food_attributes,
            energy_decrease, time_step, food_calories, reproduction_energy,
            width, height, eating_time, bird_attributes_list)

    root.mainloop()

    functions_D_A.simulation_visualization(bird_attributes_list)
