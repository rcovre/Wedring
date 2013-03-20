"""The **Wedring** error module.
This module contains the definition of the :class:`WedringError` class.
Inside the **Wedring** pipeline this class will always be raised with two
arguments: a message and a error code, respectively. **Wedring** has its own
error codes to differ from the system codes.
The error codes and their meaning are:

+----------------------------------+------+
| Status                           | Code |
+==================================+======+
| Argument list not supported      |  132 |
+----------------------------------+------+
| Insufficient number of arguments |  133 |
+----------------------------------+------+
| Invalid option value             |  134 |
+----------------------------------+------+
| Required option                  |  135 |
+----------------------------------+------+
| Invalid option                   |  136 |
+----------------------------------+------+
| File not found                   |  137 |
+----------------------------------+------+
| Program not found                |  138 |
+----------------------------------+------+
| Configuration file error         |  139 |
+----------------------------------+------+
| Unbalanced options               |  140 |
+----------------------------------+------+
| Pipeline execution error         |  141 |
+----------------------------------+------+

.. moduleauthor:: Rafael Covre <covrebio06@gmail.com>

"""

class WedringError(Exception):
    """Error class raised when something wrong happens during the **Wedring**
    execution. System errors are also raised as **Wedring** errors, using their
    error codes and messages.

    """

    def __init__(self, errno=0, strerror=""):
        """:class:WedringError constructor.

        :param message: String containing the error message.
        :type message: str
        :param code: The error code.
        :type code: int

        """
        self.errno = errno
        self.strerror = strerror
        super(WedringError, self).__init__(self.errno, self.strerror)

    def __str__(self):
        return '[Errno %d] %s.' % (self.errno, self.strerror)
