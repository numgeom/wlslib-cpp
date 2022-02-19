CXX = g++
CPPFLAGS = -DDEBUG -I.. -I../coder
CXXFLAGS = -g -rdynamic -Wall -Wextra -std=c++11
LDFLAGS =
GTEST_LIB_ROOT =
GTEST_INCLUDE_ROOT =
ifneq ($(GTEST_INCLUDE_ROOT),)
	CPPFLAGS += -I$(GTEST_INCLUDE_ROOT)
endif
ifneq ($(SANITIZE),)
	CXXFLAGS += -fno-omit-frame-pointer -fsanitize=address
endif

TEST_SRCS := $(shell \ls test_*.cpp)
TESTS = $(TEST_SRCS:%.cpp=%)
EXECS = $(TEST_SRCS:%.cpp=%.exe)
DEPFILES := $(TEST_SRCS:%.cpp=%.d)

VPATH = ../
LAPACK_LIB =
ifeq ($(LAPACK_LIB),)
	LAPACK_LIB = -llapack -lblas
endif
LIBS = $(LAPACK_LIB)
ifneq ($(GTEST_LIB_ROOT),)
	LIBS += -L$(GTEST_LIB_ROOT)
endif
LIBS += -lgtest -lgtest_main -pthread

.PHONY: clean
.PRECIOUS: $(EXECS)

all: $(TESTS)

test_%: test_%.exe
	./$<

define MAKE_EXE_RULE =
$(1).d: $(1).cpp | $(DEPDIR)
	@$(CXX) $(CPPFLAGS) -MM -MP -MT $$@ -MF $$@ $$<
include $(1).d

$(1).exe: $(1).cpp $(1).d wls.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $$< -o $$@ $(LIBS) $(LDFLAGS)
endef

$(foreach test,$(TESTS),$(eval $(call MAKE_EXE_RULE, $(test))))

clean:
	@rm -f $(EXECS) *.dSYM $(DEPFILES)

