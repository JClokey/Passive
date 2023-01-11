Usage
=====

.. _installation:

Installation
------------

To use Passive, first install it using pip:

.. code-block:: console

   (.venv) $ pip install Passive

## calculators
----------------

To calculate the sampling rate in kinetic mode,
you can use the ``passive.water.kinetic_sampling_rate()`` function:

.. autofunction:: passive.water.kinetic_sampling_rate()

The ``kind`` parameter should be either ``"meat"``, ``"fish"``,
or ``"veggies"``. Otherwise, :py:func:`lumache.get_random_ingredients`
will raise an exception.

.. autoexception:: lumache.InvalidKindError

For example:

>>> import lumache
>>> lumache.get_random_ingredients()
['shells', 'gorgonzola', 'parsley']


